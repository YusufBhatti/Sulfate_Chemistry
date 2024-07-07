! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reads in auxillary data into dump

MODULE Rcf_Aux_File_Mod

!  Subroutine Rcf_Aux_File  - reads auxillary data
!
! Description:
!   Reads data from external dumps and incorporates it into the
!   output dump. Four modes are supported
!                Tracer Data
!                User Prognostics
!                Area Tranplants
!
! Method:
!   A fields array is set up for the auxillary file and the relevant
!   fields located therein. Data is copied as appropriate for the
!   Mode (ie transplants within an area, all of a user prognostic
!   copied, upper levels only of Tracers copied).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

! Parameters describing actions possible
INTEGER, PARAMETER       :: tracers   = 1
INTEGER, PARAMETER       :: user_prog = 2
INTEGER, PARAMETER       :: transplant= 3

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_AUX_FILE_MOD'

CONTAINS

SUBROUTINE Rcf_Aux_File( Hdr_Aux, Hdr_Out, Fields_Out, Field_Count_Out,&
                         ACTION, Sctn_code_aux, Item_Code_aux,         &
                         Sctn_code_out, Item_code_out )

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE mask_compression, ONLY: compress_to_mask

USE Rcf_UMhead_Mod, ONLY: &
    Um_Header_Type

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE UM_ParCore, ONLY: &
    mype

USE rcf_trans_mod              ! all of it

USE Rcf_Write_Field_Mod, ONLY: &
    Rcf_Write_Field

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Setup_Field_mod, ONLY: &
    Rcf_Setup_Field

USE Rcf_Alloc_Field_mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Locate_mod, ONLY: &
    Rcf_Locate

USE Rcf_ReadUMhdr_Mod, ONLY: &
    Rcf_ReadUMhdr

USE Rcf_HeadAddress_Mod, ONLY:  &
IC_XLen,                         IC_YLen,    &
FH_DTYear,                       FH_VTYear,  &
FH_DTDayNo,                      FH_VTDayNo,       &
FH_Dataset,                      FH_Dataset_Ancil

USE Rcf_Grid_Type_Mod, ONLY:&
    Output_Grid

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Lsm_Mod, ONLY: &
    local_land_out,     &
    local_lsm_out

USE Rcf_Global_To_Local_Mod, ONLY: &
    Rcf_Global_To_Local_Subdomain

USE cppxref_mod, ONLY:             &
    ppx_type_real,                  &
    ppx_type_int,                   &
    ppx_type_log,                   &
    ppx_atm_compressed

USE mpp_conf_mod, ONLY: exclude_halos_ew, exclude_halos_ns

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( Um_Header_type ), INTENT(INOUT)  :: Hdr_Aux ! only unit no.,
                                                  ! file already open
TYPE( Um_Header_type ), INTENT(IN)     :: Hdr_Out
TYPE( Field_Type ), POINTER            :: Fields_Out(:)

INTEGER, INTENT(IN)                    :: Field_Count_Out
INTEGER, INTENT(IN)                    :: sctn_code_aux
INTEGER, INTENT(IN)                    :: item_code_Aux
INTEGER, INTENT(IN)                    :: sctn_code_out
INTEGER, INTENT(IN)                    :: item_code_out ! only for UP
INTEGER, INTENT(IN)                    :: ACTION

! Local variables
INTEGER                     :: i, j    ! loopers
INTEGER                     :: k,l     ! loopers
INTEGER                     :: pos_out ! field position
INTEGER                     :: pos_aux ! field position
INTEGER                     :: lrow1   ! Local positions of
INTEGER                     :: lrow2   ! trans data validity
INTEGER                     :: lcol1
INTEGER                     :: lcol2
INTEGER                     :: level_base
INTEGER                     :: level_top
INTEGER                     :: SIZE
INTEGER                     :: field_count_aux
INTEGER                     :: ErrorStatus
INTEGER                     :: Copy_Count ! counter for no. of x a field
                                          ! is copied

INTEGER, PARAMETER          :: st_no_data = -3

CHARACTER (LEN=20)          :: title
CHARACTER (LEN=*), PARAMETER:: RoutineName='RCF_AUX_FILE'
CHARACTER (LEN=errormessagelength)   :: Cmessage

TYPE( field_type ), POINTER :: fields_aux(:)

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY( fields_aux )
!------------------------------------------------------------------
! Read in the Auxillary file header
!------------------------------------------------------------------
CALL Rcf_ReadUMhdr( Hdr_Aux )

!------------------------------------------------------------------
! If tracers check that data time of the aux file matches
! verification time of the output file
!------------------------------------------------------------------
IF ( ACTION == tracers ) THEN
  DO i = 0, 6
    IF ( Hdr_Aux % FixHd( FH_DTYear + i ) /= &
         Hdr_Out % FixHd( FH_VTYear + i ) ) THEN
      WRITE(umMessage,*) 'Mismatch in date info for auxillary file'
      CALL umPrint(umMessage,src='rcf_aux_file_mod')
      WRITE(umMessage,*) 'Aux file Data times = ', (Hdr_Aux % Fixhd( j ), &
                                             j = FH_DTYear, FH_DTDayNo )
      CALL umPrint(umMessage,src='rcf_aux_file_mod')
      WRITE(umMessage,*) 'Out file Ver. times = ', (Hdr_Out % Fixhd( j ), &
                                             j = FH_VTYear, FH_VTDayNo )
      CALL umPrint(umMessage,src='rcf_aux_file_mod')

      Cmessage = 'Date information mismatch between aux and dump files'
      ErrorStatus = 10
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF
  END DO
END IF

!-------------------------------------------------------------------
! Check that resolutions of two files match
!-------------------------------------------------------------------
IF ( Hdr_Aux % IntC( IC_XLen ) /= Hdr_Out % IntC( IC_XLen ) .OR. &
     Hdr_Aux % IntC( IC_YLen ) /= Hdr_Out % IntC( IC_YLen ) ) THEN

  Cmessage = 'Dimensions of AUX file and dump file do not match'
  ErrorStatus = 20
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!-------------------------------------------------------------------
! Set up field data-types for aux file
! (Grid resolutions should be as for output grid)
!-------------------------------------------------------------------
title = 'Auxillary File'
IF ( ACTION == user_prog .AND.                                 &
     Hdr_Aux % FixHd(FH_Dataset)  == FH_Dataset_Ancil ) THEN
  ! User prognostic ancillary files have full fields for land-only
  ! fields
  CALL Rcf_Setup_Field( fields_aux, Hdr_Aux, Output_Grid,      &
                        field_count_aux, title,                &
                        Output_Grid % loc_p_rows *             &
                        Output_Grid % loc_p_row_length )
ELSE
  CALL Rcf_Setup_Field( fields_aux, Hdr_Aux, Output_Grid,      &
                        field_count_aux, title, local_land_out )
END IF

!-------------------------------------------------------------------
! Main data handling/replacement - start with TRANS data
!-------------------------------------------------------------------
IF ( ACTION == transplant ) THEN
  DO i = 1, num_trans

    IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
      WRITE(umMessage,*) 'Transplanting data for stashcode ', itemc_array( i )
      CALL umPrint(umMessage,src='rcf_aux_file_mod')
    END IF

    ! Read aux field
    CALL Rcf_Locate( sctnc_array( i ), itemc_array( i ),             &
                     fields_aux, field_count_aux, pos_aux)
    CALL Rcf_Alloc_Field( fields_aux( pos_aux ) )
    CALL Rcf_Read_Field( fields_aux( pos_aux ), Hdr_Aux,             &
                         decomp_rcf_output )

    ! Read dump field
    CALL Rcf_Locate( sctnc_array( i ), itemc_array( i ),             &
                     fields_out, field_count_out, pos_out)

    CALL Rcf_Alloc_Field( fields_out( pos_out ) )
    CALL Rcf_Read_Field( fields_out( pos_out ), Hdr_Out,             &
                         decomp_rcf_output )

    ! Cannot (yet) do a transplant for a land compressed field
    IF (fields_out(pos_out) % stashmaster % grid_type ==             &
                                            ppx_atm_compressed) THEN
      Cmessage = 'Cannot transplant a land compressed field'
      ErrorStatus = 30
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

    ! Convert the co-ordinates for replacement with local ones
    CALL Rcf_Global_to_Local_Subdomain( exclude_halos_ew, exclude_halos_ns,&
                    fields_out( pos_out ) % stashmaster % grid_type , &
                    mype,                                             &
                    Row2_array( i ), Col2_array( i ),                 &
                    Row1_array( i ), Col1_array( i ),                 &
                    lrow2, lcol2, lrow1, lcol1 )

    ! If the local area is on my pe, do the transplant
    IF ( lrow1 /= st_no_data .AND. lrow2 /= st_no_data .AND. &
         lcol1 /= st_no_data .AND. lcol2 /= st_no_data ) THEN

      SELECT CASE ( fields_out( pos_out ) % stashmaster % data_type )
      CASE (ppx_type_real)
        DO j = lrow1, lrow2
          DO k = lcol1, lcol2
            l = (j-1) * fields_out(pos_out) % row_len + k
            fields_out(pos_out) %                      &
                  DATA(l,lev1_array(i):lev2_array(i))= &
            fields_aux(pos_aux) % DATA(l,lev1_array(i):lev2_array(i))
          END DO
        END DO

      CASE (ppx_type_int)
        DO j = lrow1, lrow2
          DO k = lcol1, lcol2
            l = (j-1) * fields_out(pos_out) % row_len + k
            fields_out(pos_out) %                          &
                  data_int(l,lev1_array(i):lev2_array(i))= &
            fields_aux(pos_aux) %                          &
                  data_int(l,lev1_array(i):lev2_array(i))
          END DO
        END DO

      CASE (ppx_type_log)
        DO j = lrow1, lrow2
          DO k = lcol1, lcol2
            l = (j-1) * fields_out(pos_out) % row_len + k
            fields_out(pos_out) %                          &
                  data_log(l,lev1_array(i):lev2_array(i))= &
            fields_aux(pos_aux) %                          &
                  data_log(l,lev1_array(i):lev2_array(i))
          END DO
        END DO

      END SELECT
    END IF

    ! Write out the field
    CALL Rcf_Write_Field( fields_out( pos_out ), Hdr_Out, &
                         decomp_rcf_output)

    CALL Rcf_Dealloc_Field( fields_out( pos_out ) )
    CALL Rcf_Dealloc_Field( fields_aux( pos_aux ) )

  END DO

ELSE       ! not transplant....
  !----------------------------------------------------------------
  ! Loop through the auxillary fields for consideration for
  ! output dump inclusion.
  ! After problems with an ocean dump containing a time series diagnostic
  ! feild with identical section and item no. but different size the
  ! Copy_Count now forces only the first instance in the dump to be used
  !----------------------------------------------------------------
  Copy_Count=0
  DO i = 1, field_count_aux
    IF ( (ACTION == tracers     .OR.   &
          ACTION == user_prog ) .AND.  &
         ( item_code_aux == fields_aux( i ) % stashmaster % item .AND. &
          sctn_code_aux == fields_aux( i ) % stashmaster % section ) ) THEN

      ! Increment the copy count and check its the 1st instance.
      Copy_Count=Copy_Count + 1
      IF ( Copy_Count == 1 ) THEN

        IF ( ACTION == user_prog ) THEN
          CALL Rcf_Locate( sctn_code_out, item_code_out,               &
                           fields_out, field_count_out, pos_out)
        ELSE
          CALL Rcf_Locate( fields_aux( i ) % stashmaster % section,    &
                           fields_aux( i ) % stashmaster % item,       &
                           fields_out, field_count_out, pos_out )
        END IF

        CALL Rcf_Alloc_Field( fields_aux( i ) )
        CALL Rcf_Read_Field( fields_aux(i), Hdr_Aux, decomp_rcf_output )

        CALL Rcf_Alloc_Field( fields_out( pos_out ) )

        !---------------------------------------------------------------
        ! Copy fields if user prog
        !---------------------------------------------------------------
        IF ( ACTION == user_prog ) THEN

          ! If a land only field from ancillary, can compress onto
          ! output field. This assumes a real field only
          IF (fields_aux( i ) % stashmaster % grid_type ==             &
                                ppx_atm_compressed .AND.               &
              Hdr_Aux % FixHd( FH_Dataset) == FH_Dataset_Ancil ) THEN

            DO j = 1, fields_out( pos_out ) % levels
              CALL compress_to_mask( fields_aux( i ) % DATA(:,j),        &
                  fields_out( pos_out ) % DATA(:,j),  &
                  local_lsm_out,                      &
                  fields_aux( i ) % level_size,       &
                  SIZE )
            END DO

          ELSE    ! not land compressed

            SELECT CASE ( fields_out( pos_out ) % stashmaster % data_type)
            CASE (ppx_type_real)
              fields_out( pos_out ) % DATA( :, : ) =                   &
                                      fields_aux( i ) % DATA( :, : )

            CASE (ppx_type_int)
              fields_out( pos_out ) % Data_int( :, : ) =               &
                                     fields_aux( i ) % Data_int( :, : )

            CASE (ppx_type_log)
              fields_out( pos_out ) % Data_Log( :, : ) =               &
                                     fields_aux( i ) % Data_Log( :, : )

            END SELECT
          END IF

        ELSE             ! Must be tracers

          ! If levels don't match for tracers, issue a warning
          IF ( ACTION == tracers .AND. fields_out( pos_out ) % levels/=&
                                       fields_aux( i ) % levels) THEN
            ErrorStatus = -40
            Cmessage = 'Not all tracer levels have been initialised!'
            CALL Ereport( RoutineName, ErrorStatus, Cmessage )
          END IF

          CALL Rcf_Read_Field( fields_out( pos_out ), Hdr_Out,         &
                                                     decomp_rcf_output)

          ! Only copy the aux levels over the top of the output levels
          level_top  = fields_out( pos_out ) % levels
          level_base = fields_out( pos_out ) % levels -                &
                       fields_aux( i ) % levels + 1

          SELECT CASE ( fields_out( pos_out ) % stashmaster % data_type)
          CASE (ppx_type_real)
            fields_out( pos_out ) % DATA( :, level_base : level_top) = &
                                      fields_aux( i ) % DATA( :, : )

          CASE (ppx_type_int)
            fields_out(pos_out) % Data_int( :,level_base:level_top)= &
                                  fields_aux( i ) % Data_int( :, : )

          CASE (ppx_type_log)
            fields_out(pos_out) % Data_Log( :,level_base:level_top)= &
                                  fields_aux( i ) % Data_Log( :, : )

          END SELECT

        END IF

        CALL Rcf_Write_Field( fields_out( pos_out ), Hdr_Out,          &
                              decomp_rcf_output )

        CALL Rcf_Dealloc_Field( fields_out( pos_out ) )
        CALL Rcf_Dealloc_Field( fields_aux( i ) )

      ELSE  ! Copy_Count /= 1, must have already copied field over

        ErrorStatus = -50
        Cmessage="Was about to overwrite user_prog data with 2nd field."
        CALL EReport ( RoutineName, ErrorStatus, Cmessage )

      END IF ! Copy_Count = 1
    END IF
  END DO

END IF

! Clean up the fields
DEALLOCATE( fields_aux )
NULLIFY( fields_aux )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Aux_File
END MODULE Rcf_Aux_File_Mod
