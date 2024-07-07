! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE dust_field_conversion_mod

USE field_mod,              ONLY: field_type, ASSIGNMENT(=)
USE dustbin_conversion_mod, ONLY: convert_dust_six_to_two,                  &
                                  convert_dust_two_to_six
USE um_stashcode_mod,       ONLY: stashcode_dust1_mmr, stashcode_dust2_mmr, &
                                  stashcode_dust3_mmr, stashcode_dust4_mmr, &
                                  stashcode_dust5_mmr, stashcode_dust6_mmr
USE umPrintMgr,             ONLY: umPrint, umMessage

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! Description:
!   Subroutines to convert dust between two- and six-bin schemes and 
!   interpolate the result onto the LBC grid, then add to the output file.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DUST_FIELD_CONVERSION_MOD'

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE process_dust(j, lbc_output_control, input_file, output_file,   &
                        l_vert_interp,  interp_weights, orography_lbc_fields)

! Control routine for dust processing where bin-conversion is required.

USE fieldsfile_mod,         ONLY: fieldsfile_type
USE datafile_mod,           ONLY: datafile_type
USE lbc_output_control_mod, ONLY: lbc_output_control_type
USE interp_control_mod,     ONLY: interp_control
USE interp_weights_mod,     ONLY: interp_weights_type
USE find_fields_to_interp_mod, ONLY: fields_to_interp

IMPLICIT NONE
  
CLASS(fieldsfile_type) :: input_file
CLASS(datafile_type) :: output_file
TYPE(lbc_output_control_type) :: lbc_output_control
INTEGER :: j
LOGICAL :: l_vert_interp
TYPE(interp_weights_type) :: interp_weights(:,:)
TYPE(field_type) :: orography_lbc_fields(:,:)
  
INTEGER :: new_field_num, p
TYPE(field_type), ALLOCATABLE :: input_dust_fields(:), output_dust_fields(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'PROCESS_DUST'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate the fields for both input and output dust fields
ALLOCATE(input_dust_fields(input_file%number_of_dust_bins_in_file()))
ALLOCATE(output_dust_fields(lbc_output_control%num_dust_fields))

! Assign the input dust fields
CALL input_file%read_field(fields_to_interp(j))
input_dust_fields(1) = input_file%fields(fields_to_interp(j))
CALL input_file%unload_field(fields_to_interp(j))

CALL input_file%read_field(fields_to_interp(j+1))
input_dust_fields(2) = input_file%fields(fields_to_interp(j+1))
CALL input_file%unload_field(fields_to_interp(j+1))

IF (input_file%number_of_dust_bins_in_file() == 6) THEN

  CALL input_file%read_field(fields_to_interp(j+2))
  input_dust_fields(3) = input_file%fields(fields_to_interp(j+2))
  CALL input_file%unload_field(fields_to_interp(j+2))

  CALL input_file%read_field(fields_to_interp(j+3))
  input_dust_fields(4) = input_file%fields(fields_to_interp(j+3))
  CALL input_file%unload_field(fields_to_interp(j+3))

  CALL input_file%read_field(fields_to_interp(j+4))
  input_dust_fields(5) = input_file%fields(fields_to_interp(j+4))
  CALL input_file%unload_field(fields_to_interp(j+4))

  CALL input_file%read_field(fields_to_interp(j+5))
  input_dust_fields(6) = input_file%fields(fields_to_interp(j+5))
  CALL input_file%unload_field(fields_to_interp(j+5))

END IF

WRITE(umMessage,'(A,I1,A,I1,A,A)') '[INFO] Converting ',                             &
               SIZE(input_dust_fields), ' bin dust to ',                             &
               SIZE(output_dust_fields), ' bin dust for time ',                      &
               input_dust_fields(1)%return_validity_time_string()
CALL umPrint(umMessage, src='process_dust')

! Convert between dust bins
CALL convert_dust_fields(input_dust_fields, output_dust_fields)

! Interpolate the new bins and add to the output file
DO p = 1, SIZE(output_dust_fields)
  new_field_num = output_file%add_field(interp_control(lbc_output_control, &
                                    input_file, interp_weights,                      &
                                    orography_lbc_fields, l_vert_interp,             &
                                    field_in=output_dust_fields(p)))
  CALL output_dust_fields(p)%unload_data()
  WRITE(umMessage,'(A,I5,A,A)') '[INFO] Generating '                                 &
          //'field for STASH ',                                                      & 
          output_file%fields(new_field_num)%quantity_ident, ' at time ',             &
          output_file%fields(new_field_num)%return_validity_time_string()
  CALL umPrint(umMessage, src='process_dust')

  CALL output_file%write_field(new_field_num)

  CALL output_file%fields(new_field_num)%unload_data()

END DO ! Loop over dust in output file

DEALLOCATE(output_dust_fields)
DEALLOCATE(input_dust_fields)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE process_dust

!-------------------------------------------------------------------------------

SUBROUTINE convert_dust_fields(input_fields, output_fields)
! Convert dust fields between two-bin and six-bin schemes. Wraps a UM routine
! which requires 1D input arrays; we use Field objects in CreateBC so this
! requires a format conversion.
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
TYPE(field_type) :: input_fields(:), output_fields(:)

INTEGER :: ndust_in, ndust_out
INTEGER :: data_size, nx, ny, nz
REAL, ALLOCATABLE ::    dust_div1(:), dust_div2(:), dust_div3(:),          &
            dust_div4(:), dust_div5(:), dust_div6(:)

! Error reporting
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'CONVERT_DUST_FIELDS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ndust_in = SIZE(input_fields)
ndust_out = SIZE(output_fields)

nx = SIZE(input_fields(1)%rdata(:,:,:),1)
ny = SIZE(input_fields(1)%rdata(:,:,:),2)
nz = SIZE(input_fields(1)%rdata(:,:,:),3)
data_size = nx * ny * nz

IF (ndust_in == ndust_out) THEN
  output_fields = input_fields
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
  
ALLOCATE(dust_div1(data_size))
ALLOCATE(dust_div2(data_size))
ALLOCATE(dust_div3(data_size))
ALLOCATE(dust_div4(data_size))
ALLOCATE(dust_div5(data_size))
ALLOCATE(dust_div6(data_size))

IF (ndust_in == 6 .AND. ndust_out == 2) THEN
    
  CALL convert_rdata_to_1d(input_fields(1), dust_div1)
  CALL convert_rdata_to_1d(input_fields(2), dust_div2)
  CALL convert_rdata_to_1d(input_fields(3), dust_div3)
  CALL convert_rdata_to_1d(input_fields(4), dust_div4)
  CALL convert_rdata_to_1d(input_fields(5), dust_div5)
  CALL convert_rdata_to_1d(input_fields(6), dust_div6)

  output_fields(1) = input_fields(1)
  output_fields(2) = input_fields(2)

  CALL input_fields(1)%unload_data()
  CALL input_fields(2)%unload_data()
  CALL input_fields(3)%unload_data()
  CALL input_fields(4)%unload_data()
  CALL input_fields(5)%unload_data()
  CALL input_fields(6)%unload_data()

  CALL convert_dust_six_to_two(data_size, dust_div1, dust_div2, dust_div3, &
                               dust_div4, dust_div5, dust_div6)

  CALL convert_1d_to_rdata(dust_div1, output_fields(1))
  CALL convert_1d_to_rdata(dust_div2, output_fields(2))

  output_fields(1)%quantity_ident = stashcode_dust1_mmr
  output_fields(2)%quantity_ident = stashcode_dust2_mmr

ELSE IF (ndust_in == 2 .AND. ndust_out == 6) THEN
  CALL convert_rdata_to_1d(input_fields(1), dust_div1)
  CALL convert_rdata_to_1d(input_fields(2), dust_div2)

  output_fields(1) = input_fields(1)
  output_fields(2) = input_fields(1)
  output_fields(3) = input_fields(1)
  output_fields(4) = input_fields(1)
  output_fields(5) = input_fields(1)
  output_fields(6) = input_fields(1)

  CALL input_fields(1)%unload_data()
  CALL input_fields(2)%unload_data()

  CALL convert_dust_two_to_six(data_size, dust_div1, dust_div2, dust_div3, &
                               dust_div4, dust_div5, dust_div6)

  CALL convert_1d_to_rdata(dust_div1, output_fields(1))
  CALL convert_1d_to_rdata(dust_div2, output_fields(2))
  CALL convert_1d_to_rdata(dust_div3, output_fields(3))
  CALL convert_1d_to_rdata(dust_div4, output_fields(4))
  CALL convert_1d_to_rdata(dust_div5, output_fields(5))
  CALL convert_1d_to_rdata(dust_div6, output_fields(6))

  output_fields(1)%quantity_ident = stashcode_dust1_mmr
  output_fields(2)%quantity_ident = stashcode_dust2_mmr
  output_fields(3)%quantity_ident = stashcode_dust3_mmr
  output_fields(4)%quantity_ident = stashcode_dust4_mmr
  output_fields(5)%quantity_ident = stashcode_dust5_mmr
  output_fields(6)%quantity_ident = stashcode_dust6_mmr

ELSE
  WRITE(cmessage,'(A,I1,A,I1,A)') 'Invalid dust bin conversion: ',ndust_in,&
              ' bin to ',ndust_out,' bin not supported'
  icode = 15
  CALL ereport(routinename, icode, cmessage)  

END IF

DEALLOCATE(dust_div6)
DEALLOCATE(dust_div5)
DEALLOCATE(dust_div4)
DEALLOCATE(dust_div3)
DEALLOCATE(dust_div2)
DEALLOCATE(dust_div1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE convert_dust_fields

!-------------------------------------------------------------------------------

SUBROUTINE convert_rdata_to_1d(field, output_array)
! Convert the real data in a field to a 1D array
IMPLICIT NONE
TYPE(field_type) :: field
REAL :: output_array(:)

INTEGER :: nx, ny, nz
INTEGER :: i,j,k,l
CHARACTER(LEN=*), PARAMETER :: routinename = 'CONVERT_RDATA_TO_1D'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nx = SIZE(field%rdata(:,:,:),1)
ny = SIZE(field%rdata(:,:,:),2)
nz = SIZE(field%rdata(:,:,:),3)

l = 0
DO k = 1, nz
  DO j = 1, ny
    DO i = 1, nx
      l = l + 1
      output_array(l) = field%rdata(i,j,k)
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE convert_rdata_to_1d

!-------------------------------------------------------------------------------

SUBROUTINE convert_1d_to_rdata(input_array, field)
! Copy a 1D data array into the real data storage of a field object
IMPLICIT NONE
TYPE(field_type) :: field
REAL :: input_array(:)

INTEGER :: nx, ny, nz
INTEGER :: i,j,k,l
CHARACTER(LEN=*), PARAMETER :: routinename = 'CONVERT_1D_TO_RDATA'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

nx = SIZE(field%rdata(:,:,:),1)
ny = SIZE(field%rdata(:,:,:),2)
nz = SIZE(field%rdata(:,:,:),3)

l = 0
DO k = 1, nz
  DO j = 1, ny
    DO i = 1, nx
      l = l + 1
      field%rdata(i,j,k) = input_array(l)
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE convert_1d_to_rdata

!-------------------------------------------------------------------------------

END MODULE dust_field_conversion_mod

