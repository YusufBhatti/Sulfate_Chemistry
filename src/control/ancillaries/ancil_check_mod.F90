! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Module: Contains routines used to check validity of ancils
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Ancillaries
!
MODULE ancil_check_mod

USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim
USE ereport_mod,  ONLY: ereport
USE umPrintMgr,   ONLY: umPrint, umMessage, newline
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

PRIVATE

PUBLIC :: ancil_check_grid_stagger, ancil_check_horizontal_grid

CHARACTER (LEN=*), PARAMETER :: ModuleName = 'ANCIL_CHECK_MOD'
! DrHook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CONTAINS

!--------------------------------------------------------------------------

SUBROUTINE ancil_check_horizontal_grid(lookup,                        &
                                       ancil_file,    len1_lookup,    &
                                       p_rows,        p_row_length,   &
                                       u_rows,        u_row_length,   &
                                       v_rows,        v_row_length,   &
                                       river_rows,    river_row_length)
USE ancil_mod, ONLY: ancil_file_type
USE lookup_addresses, ONLY: item_code, lbrow, lbnpt, lbpack
USE cppxref_mod, ONLY:                               &
    ppx_grid_type,                                   &    
   ! Indicators for horiz grid type
    u_points => ppx_atm_cuall,                       &
    v_points => ppx_atm_cvall,                       &
    p_points => ppx_atm_tall,                        & 
    ozone_points =>  ppx_atm_ozone,                  &
    river_routing_points => ppx_atm_river,           &
    compressed_land_points => ppx_atm_compressed,    &
    p_points_values_on_sea_only => ppx_atm_tsea
USE um_parcore, ONLY: mype
USE ozone_inputs_mod, ONLY: zon_av_ozone, zon_av_tpps_ozone

USE Packing_Codes_Mod, ONLY:PC_No_CompressType

!
! Description: 
!   Checks the ancilary file horizontal grid lookup values conform to the model 
!   grid.
! Method:
!   Takes the lookup array of the current ancillary and loops over all fields
!   If the field in the lookup matches the requested ancil stashcodes from
!   the current file then grid checks will be applied.
IMPLICIT NONE
     
! Arguments
INTEGER, INTENT(IN) :: lookup(:,:)
TYPE (ancil_file_type), INTENT(IN) :: ancil_file ! Ancil file type
INTEGER, INTENT(IN) :: len1_lookup
INTEGER, INTENT(IN) :: p_rows       ! No. of rows for pressure-type variables
INTEGER, INTENT(IN) :: p_row_length ! Row length for pressure-type variables   
INTEGER, INTENT(IN) :: u_rows       ! No. of rows for u wind-type variables 
INTEGER, INTENT(IN) :: u_row_length ! Row length for u wind type variables  
INTEGER, INTENT(IN) :: v_rows       ! No. of rows for v wind-type variables 
INTEGER, INTENT(IN) :: v_row_length ! Row length for v wind-type variables  
INTEGER, INTENT(IN) :: river_rows   ! No. of rows for river routing-type 
                                    ! variables
INTEGER, INTENT(IN) :: river_row_length ! Row length for river routing-type 
                                        ! variables  

! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'ANCIL_CHECK_HORIZONTAL_GRID'
CHARACTER (LEN=errormessagelength*2) :: cmessageLong  ! Error message
INTEGER :: icode
INTEGER :: i_stash   ! loop variable for stashcodes in file
INTEGER :: stashcode
INTEGER :: i_lookup
INTEGER :: num_lookups
INTEGER :: grid_type_code

! DrHook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

num_lookups = SIZE(lookup,2)
! Loop over all stash that have been requested 
DO i_stash = 1, ancil_file%num_stash
  stashcode = ancil_file%stashcodes(i_stash)
  grid_type_code = get_stashmaster_element(stashcode, ppx_grid_type)

  ! Loop over all lookup fields in file and look for match
  DO i_lookup = 1, num_lookups
    IF (stashcode == lookup(item_code, i_lookup)) THEN

      SELECT CASE(grid_type_code)

      CASE(p_points, p_points_values_on_sea_only)
        CALL check_lookup_values_int(lookup(lbrow, i_lookup), p_rows,         &
             "p rows", "")
        CALL check_lookup_values_int(lookup(lbnpt, i_lookup), p_row_length,   &
             "p row length", "")
        
      CASE(u_points)
        CALL check_lookup_values_int(lookup(lbrow, i_lookup), u_rows,         &
              "u rows", "")
        CALL check_lookup_values_int(lookup(lbnpt, i_lookup), u_row_length,   &
              "u row length", "")

      CASE(v_points)
        CALL check_lookup_values_int(lookup(lbrow, i_lookup), v_rows,         &
              "v rows", "")
        CALL check_lookup_values_int(lookup(lbnpt, i_lookup), v_row_length,   &
              "v row length", "")
        
      CASE(ozone_points)
        ! Check whether zonal or full field ozone as this will effect the 
        ! expected row length
        IF (zon_av_ozone .OR. zon_av_tpps_ozone) THEN
          CALL check_lookup_values_int(lookup(lbnpt, i_lookup), 1,            &
               "ozone row length",   newline //                               &
               "ZONAL ozone selected, row length should equal 1")
        ELSE
          CALL check_lookup_values_int(lookup(lbnpt, i_lookup), p_row_length, &
               "ozone row length",   newline //                               &
               "Full field ozone selected, row length should match " //       &
               "p row length")
        END IF
        
        CALL check_lookup_values_int(lookup(lbrow, i_lookup), p_rows,         &
             "ozone rows", newline // "Ozone rows should match p rows")

      CASE(river_routing_points)
        CALL check_lookup_values_int(lookup(lbrow, i_lookup), river_rows,     &
              "river rows", "")
        CALL check_lookup_values_int(lookup(lbnpt, i_lookup),                 &
             river_row_length, "river row length", "")

      CASE(compressed_land_points)
        IF ( MOD( lookup(lbpack, i_lookup)/10, 10 ) == PC_No_CompressType ) THEN
          ! If the STASHmaster specifies that the field is on compressed land 
          ! point but the field is not packed then the UM allows the field to 
          ! be on the standard non compressed grid
          CALL check_lookup_values_int(lookup(lbrow, i_lookup), p_rows,       &
              "p rows",  newline // "Field grid type is compressed " //       &
              "land points, but data is uncompressed." // newline //          &
              "In this case the field should match the p grid."  )
          CALL check_lookup_values_int(lookup(lbnpt, i_lookup), p_row_length, &
              "p row length", newline // "Field grid type is compressed "//   &
              "land points, but data is uncompressed." // newline //          &
              "In this case the field should match the p grid.")
        ELSE
          ! Otherwise land packed fields should have rows and row length zero
          CALL check_lookup_values_int(lookup(lbrow, i_lookup), 0,            &
              "rows", newline // "Land compressed fields lookup rows should " &
              // "equal 0")
          CALL check_lookup_values_int(lookup(lbnpt, i_lookup), 0,            &
              "row length", newline // "Land compressed fields row "   &
              // "length equal 0")
        END IF

       CASE DEFAULT
         WRITE(cmessageLong, '(A,I0,A,I0,A,I0)')                              &
              "Unsupported grid type code : ", grid_type_code, newline //     &
              TRIM(print_current_ancil_info())
         icode = 20
         CALL ereport(routinename, icode, cmessageLong)
      END SELECT
        
    END IF ! (stashcode == lookup(item_code, i_lookup))
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

CONTAINS 

!----------------------------

SUBROUTINE check_lookup_values_int(ancil_lookup_value, expected_value, &
     name_of_check, additional_error_info)
! Internal routine used to compare two values and return an ereport if the
! values do not match. Expects the name of the check to be passed in as well
! as any additional error information that needs to be passed to the 
! error message.
IMPLICIT NONE

INTEGER, INTENT(IN) :: expected_value
INTEGER, INTENT(IN) :: ancil_lookup_value
CHARACTER (LEN=*), INTENT(IN) :: name_of_check
CHARACTER (LEN=*), INTENT(IN) :: additional_error_info


IF (expected_value /= ancil_lookup_value) THEN
  WRITE(cmessageLong, '(A,I0,A,I0,A,I0,A,I0)')                               &
       "Mismatch between model and ancil field " // name_of_check            &
       // TRIM(additional_error_info) // newline //                          &
       "Expected " // name_of_check // " : ", expected_value, newline//      &
       "Ancil    " // name_of_check // " : ", ancil_lookup_value, newline // &
       TRIM(print_current_ancil_info())
  icode = 10
  CALL ereport(routinename, icode, cmessageLong)
END IF
END SUBROUTINE check_lookup_values_int

!----------------------------

FUNCTION print_current_ancil_info()
! Internal function used to return current ancil info as a character string
! Operates directly on the variables in the containing subroutine
IMPLICIT NONE
CHARACTER (LEN=errormessagelength) :: print_current_ancil_info
WRITE(print_current_ancil_info, '(A,I0,A,I0)') &
       "Ancil file : " //  ancil_file % filename  //   newline //            &
       "Lookup num : ", i_lookup,                      newline //            &
       "Stashcode  : ", stashcode
END FUNCTION print_current_ancil_info

END SUBROUTINE ancil_check_horizontal_grid

!----------------------------------------------------------------------------

INTEGER FUNCTION get_stashmaster_element(stashcode, element)
#if defined(RECON)
USE rcf_exppx_mod, ONLY: rcf_exppx
USE rcf_ppx_info_mod, ONLY: stm_record_type
#endif

USE ppxlook_mod,        ONLY: exppxi
USE cppxref_mod,        ONLY: ppx_grid_type    
USE missing_data_mod,   ONLY: imdi    

! We need method of obtaining the stashmaster record information
! that is independent of whether we are running this code in the 
! UM or recon. Unfortunately the UM and recon have different methods
! of reading and storing the STASHmaster.  Until such a point where
! this is rationalised (if ever) it will be necessary to use this 
! routine as an interface to the UM and recon methods.
!
! Given a stashcode and element (e.g. grid code), returns the value of
! that element for that stashcode.
IMPLICIT NONE

INTEGER, INTENT(IN) :: stashcode, element

! Local variables
INTEGER :: section, item
INTEGER, PARAMETER ::  model = 1   ! Atmosphere model only
INTEGER :: icode
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'GET_STASHMASTER_ELEMENT'
CHARACTER (LEN=errormessagelength) :: cmessage  ! Error message

#if defined(RECON)
TYPE (stm_record_type), POINTER ::  stm_record
#endif 

! DrHook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

section = stashcode / 1000
item = MOD(stashcode, 1000)

!------------- RECON ONLY CODE --------------------------------------------
#if defined(RECON)
! Passing in NoFindArg=false means the recon routine will ereport with a fatal
! error if the STASHmaster record is not found.
stm_record => rcf_exppx(model, section, item, NoFindArg=.FALSE.)
IF (element == ppx_grid_type) THEN
  get_stashmaster_element = stm_record % grid_type
END IF
!------------- END OF RECON ONLY CODE -------------------------------------

!------------- UM ONLY CODE -----------------------------------------------
#else 
get_stashmaster_element = exppxi(model, section, item, element)
! UM routine will not ereport if the stashmaster record isn't found - it will
! produce an imdi value. Use that to ereport with an error here instead.
IF (get_stashmaster_element == imdi) THEN
  WRITE (cmessage, '(A, I0, A, I0, A)')                       &
       'Cant find required STASH item ', item, ' section ', section, &
       ' in STASHmaster'
  CALL ereport(routinename, icode, cmessage)
END IF
#endif
!------------- END OF UM ONLY CODE ----------------------------------------


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION get_stashmaster_element

!--------------------------------------------------------------------------

SUBROUTINE ancil_check_grid_stagger(model_grid_stagger, ancil_grid_stagger, &
                                    ancil_file)
USE ancil_mod, ONLY: ancil_file_type

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN) :: model_grid_stagger
INTEGER, INTENT(IN) :: ancil_grid_stagger
TYPE(ancil_file_type), INTENT(IN) :: ancil_file

! Local variables
INTEGER :: icode
CHARACTER(LEN=1024) :: cmessageLong  
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'ANCIL_CHECK_GRID_STAGGER'

! DrHook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (model_grid_stagger /= ancil_grid_stagger) THEN
  icode=-10
  WRITE(cmessageLong, '(A,I0,A,I0,A)')                                        & 
     "Ancil file mismatch in fixed header(9) grid stagger value " //newline// &
     "Model grid stagger = ", model_grid_stagger,                   newline// &
     "Ancil file grid stagger = ", ancil_grid_stagger,              newline// &
     "Ancil file path = " // ancil_file%filename                  //newline// &
     "PLEASE READ - this warning will be converted to an error"   //newline// &
     "in future. Please update ancil file to specify the correct" //newline// &
     "grid stagger value."
  CALL ereport(routinename, icode, cmessageLong )
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ancil_check_grid_stagger

!--------------------------------------------------------------------------

END MODULE ancil_check_mod
