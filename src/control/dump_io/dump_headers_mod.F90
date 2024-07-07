! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Contains infomation related to the dump headers
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Dump I/O

MODULE dump_headers_mod

USE parkind1,             ONLY: jprb, jpim
USE yomhook,              ONLY: lhook, dr_hook
USE ereport_mod,          ONLY: ereport
USE umPrintMgr,           ONLY: umMessage, newline, umPrint
USE nlsizes_namelist_mod, ONLY: len_fixhd, a_len_inthd, a_len_cfi1,   &
                                a_len_cfi2, a_len_cfi3, a_len_realhd, &
                                a_len1_levdepc, a_len2_levdepc,       &
                                a_len1_rowdepc, a_len2_rowdepc,       &
                                a_len1_coldepc, a_len2_coldepc,       &
                                a_len1_flddepc, a_len2_flddepc,       &
                                a_len_extcnst, len_dumphist,          &
                                len1_lookup, a_len2_lookup,           &
                                mpp_len1_lookup, a_len2_lookup
USE atm_fields_bounds_mod, ONLY: tdims
USE missing_data_mod, ONLY: rmdi

USE stochastic_physics_run_mod,  ONLY: l_skeb2, l_rp2, l_spt,         &
    stph_n1, stph_n2, stph_spt_data_check, stph_spt_data_present,     &
    stph_skeb2_data_check, stph_skeb2_data_present,                   &
    stph_rp2_data_check, stph_rp2_data_present, stph_header_flag,     &
    rp_max

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: modulename = 'DUMP_HEADERS_MOD'

! Meaningful parameter names for integer constants header
INTEGER, PARAMETER :: ih_a_step          = 1  ! Timestep no.
INTEGER, PARAMETER :: ih_rowlength       = 6  ! No. of points E-W
INTEGER, PARAMETER :: ih_rows            = 7  ! No. of points N-S

! No. of model levels (0=surface)
INTEGER, PARAMETER :: ih_model_levels    = 8

! No. of model levels with moisture
INTEGER, PARAMETER :: ih_wet_levels      = 9

! No. of deep soil temperature levels
INTEGER, PARAMETER :: ih_soilT_levels    = 10

INTEGER, PARAMETER :: ih_cloud_levels    = 11 ! No. of cloud levels
INTEGER, PARAMETER :: ih_tracer_levels   = 12 ! No. of tracer levels

! No. of boundary layer levels
INTEGER, PARAMETER :: ih_boundary_levels = 13
INTEGER, PARAMETER :: ih_N_types         = 15 ! No. of field types

! Height generation method
INTEGER, PARAMETER :: ih_height_gen      = 17


! First rho level at which height is constant
INTEGER, PARAMETER :: ih_1_c_rho_level   = 24
INTEGER, PARAMETER :: ih_land_points     = 25 ! No. of land points
INTEGER, PARAMETER :: ih_ozone_levels    = 26 ! No. of ozone levels

! No. of deep soil moisture levels
INTEGER, PARAMETER :: ih_soilQ_levels    = 28

! Stochastic physics values
INTEGER, PARAMETER :: ih_stochastic_flag = 29 ! Flag for presence of SP data
INTEGER, PARAMETER :: ih_stph_n1         = 30 ! First dimension of SP data
INTEGER, PARAMETER :: ih_stph_n2         = 31 ! Second dimension of SP data
INTEGER, PARAMETER :: ih_sp_seed         = 32 ! Seed for SP random number 
                                              ! generator

! Number of convective cloud levels
INTEGER, PARAMETER :: ih_convect_levels  = 34
INTEGER, PARAMETER :: ih_rad_step        = 35 ! Radiation timestep
INTEGER, PARAMETER :: ih_AMIP_flag       = 36 ! Flag for AMIP run
INTEGER, PARAMETER :: ih_AMIP_year       = 37 ! First AMIP year
INTEGER, PARAMETER :: ih_AMIP_month      = 38 ! First AMIP month
INTEGER, PARAMETER :: ih_AMIP_day        = 49 ! First AMIP day
INTEGER, PARAMETER :: ih_ozone_month     = 40 ! Current ozone month
INTEGER, PARAMETER :: ih_SH_zonal_quad   = 41 ! L_SH_zonal_quadratics
INTEGER, PARAMETER :: ih_SH_zonal_begin  = 42 ! SH_zonal_begin
INTEGER, PARAMETER :: ih_SH_zonal_period = 43 ! SH_zonal_period
INTEGER, PARAMETER :: ih_SH_level_weight = 44 ! SuHe_level_weight
INTEGER, PARAMETER :: ih_SH_sigma_cutoff = 45 ! SuHe_sigma_cutoff
INTEGER, PARAMETER :: ih_friction_time   = 46 ! frictional_timescale

! Meaningful parameter names for real constants header

! East-West   grid spacing in degrees
INTEGER, PARAMETER :: rh_deltaEW         = 1

! North-South grid spacing in degrees
INTEGER, PARAMETER :: rh_deltaNS         = 2

! Latitude  of first p point in degrees
INTEGER, PARAMETER :: rh_baselat         = 3

! Longitude of first p point in degrees
INTEGER, PARAMETER :: rh_baselong        = 4

! Latitude  of rotated N pole in degrees
INTEGER, PARAMETER :: rh_rotlat          = 5

! Longitude of rotated N pole in degrees
INTEGER, PARAMETER :: rh_rotlong         = 6

! Height of top theta level (m)
INTEGER, PARAMETER :: rh_z_top_theta     =16

! total moisture of the atmosphere
INTEGER, PARAMETER :: rh_tot_m_init      =18

! total mass of atmosphere
INTEGER, PARAMETER :: rh_tot_mass_init   =19

! total energy of atmosphere
INTEGER, PARAMETER :: rh_tot_energy_init =20

! energy correction = energy drift
INTEGER, PARAMETER :: rh_energy_corr     =21

! Meaningful parameter names for fixed header

! Grid co-ordinate system (e.g. Long-Lat)
INTEGER, PARAMETER :: fh_CoordSystem    = 11
! Recognised values for fh_CoordSystem (IMDI is assumed long-lat):
INTEGER, PARAMETER :: fh_CoordSystem_RegLonLat = 1
INTEGER, PARAMETER :: fh_CoordSystem_Cartesian = 9

! Start of Row Dependent Constant
INTEGER, PARAMETER :: fh_RowDepCStart   = 115

! Start of Col Dependent Constant
INTEGER, PARAMETER :: fh_ColDepCStart   = 120

!Field-dependant constants data positions
INTEGER :: fdc_skeb2_dpsidtc_start, fdc_skeb2_dpsidtc_length
INTEGER :: fdc_skeb2_dpsidts_start, fdc_skeb2_dpsidts_length
INTEGER :: fdc_spt_coeffc_start, fdc_spt_coeffc_length
INTEGER :: fdc_spt_coeffs_start, fdc_spt_coeffs_length
INTEGER :: fdc_rp2_coef_start, fdc_rp2_coef_length


! Arrays

INTEGER :: a_fixhd(len_fixhd)      ! fixed length header

INTEGER, ALLOCATABLE :: a_inthd(:) ! integer header
INTEGER, ALLOCATABLE :: a_cfi1(:)  ! compress field index
INTEGER, ALLOCATABLE :: a_cfi2(:)  ! compress field index
INTEGER, ALLOCATABLE :: a_cfi3(:)  ! compress field index

! lookup headers
INTEGER, ALLOCATABLE :: a_lookup(:,:)    
INTEGER, ALLOCATABLE :: a_mpp_lookup(:,:)

REAL, ALLOCATABLE :: a_realhd(:)     ! real header
REAL, ALLOCATABLE :: a_levdepc(:)    ! level dep const
REAL, ALLOCATABLE :: a_rowdepc(:)    ! row dep const
REAL, ALLOCATABLE :: a_coldepc(:)    ! column dep const
REAL, ALLOCATABLE :: a_flddepc_in(:) ! field dep const (in)
REAL, ALLOCATABLE :: a_flddepc(:)    ! field dep const
REAL, ALLOCATABLE :: a_extcnst(:)    ! extra constants
REAL :: a_dumphist(len_dumphist+1)   ! temp hist file
LOGICAL, ALLOCATABLE :: sea_mask(:,:)

CONTAINS

SUBROUTINE allocate_dump_headers()

IMPLICIT NONE

INTEGER                       :: icode
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: routinename = 'ALLOCATE_DUMP_HEADERS'

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_in,zhook_handle)

icode=0

IF (ALLOCATED(a_inthd     )) DEALLOCATE(a_inthd)
IF (ALLOCATED(a_cfi1      )) DEALLOCATE(a_cfi1 )
IF (ALLOCATED(a_cfi2      )) DEALLOCATE(a_cfi2 )
IF (ALLOCATED(a_cfi3      )) DEALLOCATE(a_cfi3 )
IF (ALLOCATED(a_lookup    )) DEALLOCATE(a_lookup)
IF (ALLOCATED(a_mpp_lookup)) DEALLOCATE(a_mpp_lookup)
IF (ALLOCATED(a_realhd    )) DEALLOCATE(a_realhd)
IF (ALLOCATED(a_levdepc   )) DEALLOCATE(a_levdepc)
IF (ALLOCATED(a_rowdepc   )) DEALLOCATE(a_rowdepc)
IF (ALLOCATED(a_coldepc   )) DEALLOCATE(a_coldepc)
IF (ALLOCATED(a_flddepc   )) DEALLOCATE(a_flddepc)
IF (ALLOCATED(a_flddepc_in)) DEALLOCATE(a_flddepc_in)
IF (ALLOCATED(a_extcnst   )) DEALLOCATE(a_extcnst)
IF (ALLOCATED(sea_mask    )) DEALLOCATE(sea_mask)

ALLOCATE( a_inthd(a_len_inthd), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_inthd]'//        newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_cfi1(a_len_cfi1+1), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_cfi1]'//         newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_cfi2(a_len_cfi2+1), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_cfi2]'//         newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_cfi3(a_len_cfi3+1), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_cfi3]'//         newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_lookup(len1_lookup,a_len2_lookup), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_lookup]'//       newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_mpp_lookup(mpp_len1_lookup,a_len2_lookup), STAT=icode,           &
          ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_mpp_lookup]'//   newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_realhd(a_len_realhd), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_realhd]'//       newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_levdepc(a_len1_levdepc*a_len2_levdepc+1), STAT=icode,            &
          ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_levdepc]'//      newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_rowdepc(a_len1_rowdepc*a_len2_rowdepc+1), STAT=icode,            &
          ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_rowdepc]'//      newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_coldepc(a_len1_coldepc*a_len2_coldepc+1), STAT=icode,            &
          ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_coldepc]'//      newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_flddepc(a_len1_flddepc*a_len2_flddepc+1), STAT=icode,            &
          ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_flddepc]'//   newline//    &
                TRIM(umMessage) )
END IF

ALLOCATE( a_flddepc_in(a_len1_flddepc*a_len2_flddepc+1), STAT=icode,         &
          ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_flddepc_in]'//   newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( a_extcnst(a_len_extcnst+1), STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_extcnst]'//      newline// &
                TRIM(umMessage) )
END IF

ALLOCATE( sea_mask(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end),    &
          STAT=icode, ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [sea_mask]'//       newline// &
                TRIM(umMessage) )
END IF

IF (lhook) CALL dr_hook(modulename//':'//routinename,zhook_out,zhook_handle)

RETURN

END SUBROUTINE allocate_dump_headers


SUBROUTINE allocate_sp_coefficients()

IMPLICIT NONE

INTEGER             :: fdc_size, fdc_start, spt_size, skeb2_size, rp_size
INTEGER             :: ix_i, ix_j

fdc_size = 1

! size required for sine and cosine coefficients for SPT scheme
IF ( l_spt ) THEN
  spt_size = (stph_n2 + 1) * stph_n2
  stph_header_flag = stph_header_flag + stph_spt_data_present
ELSE
  spt_size = 0
END IF

! size required for sine and cosine coefficients for SKEB2 scheme
IF ( l_skeb2 ) THEN
  skeb2_size = (stph_n2 + 1) * (stph_n2 - stph_n1 + 1)
  stph_header_flag = stph_header_flag + stph_skeb2_data_present
ELSE
  skeb2_size = 0
END IF

! Size required for the RP coefficients in the RP2 scheme
IF ( l_rp2 ) THEN
  rp_size = rp_max
  stph_header_flag = stph_header_flag + stph_rp2_data_present
ELSE
  rp_size = 0
END IF

CALL calculate_fdc_positions(spt_size, skeb2_size, rp_size, fdc_size)

! total size to allocate in field of constants
fdc_size = fdc_size +                 &
           fdc_spt_coeffc_length +    &
           fdc_spt_coeffs_length +    &
           fdc_skeb2_dpsidtc_length + &
           fdc_skeb2_dpsidts_length + &
           fdc_rp2_coef_length


a_len1_flddepc = fdc_size
a_len2_flddepc = 1

! calculate the start position of the field-dependant constants from 
! other entries
IF (a_fixhd(125) > 0) THEN
  fdc_start = a_fixhd(125)
ELSE
  fdc_start = calculate_fdc_start()
END IF

! set the start and size of the field-dependant constants in the fixed header
 a_fixhd(125) = fdc_start
 a_fixhd(126) = a_len1_flddepc
 a_fixhd(127) = a_len2_flddepc

! offset the start position of dump elements that come after field-dependant
! constants
IF (a_fixhd(140) > 0) THEN
  a_fixhd(140) = a_fixhd(140) + fdc_size ! compressed field 1 start
END IF

IF (a_fixhd(142) > 0) THEN
  a_fixhd(142) = a_fixhd(142) + fdc_size ! compressed field 2 start
END IF

IF (a_fixhd(144) > 0) THEN
  a_fixhd(144) = a_fixhd(144) + fdc_size ! compressed field 3 start
END IF

IF (a_fixhd(150) > 0) THEN
  a_fixhd(150) = a_fixhd(150) + fdc_size ! lookup start
END IF

IF (a_fixhd(160) > 0) THEN
  a_fixhd(160) = a_fixhd(160) + fdc_size ! data start
END IF


! allocate the new field-dependant constants array
IF ( ALLOCATED(a_flddepc) ) THEN
  DEALLOCATE(a_flddepc)
END IF

ALLOCATE( a_flddepc(a_len1_flddepc * a_len2_flddepc + 1) )

! initialise the values of the field-dependant constants array
DO ix_j = 1,a_len2_flddepc
  DO ix_i =  1,a_len1_flddepc
    a_flddepc( (ix_j-1) * a_len1_flddepc + ix_i ) = rmdi
  END DO
END DO

END SUBROUTINE allocate_sp_coefficients

SUBROUTINE calculate_fdc_positions(spt_size, skeb2_size, rp_size, fdc_size)

IMPLICIT NONE

INTEGER, INTENT(IN) :: spt_size, skeb2_size, rp_size, fdc_size

! calculate sizes and positions of SPT data in header array
fdc_spt_coeffc_start = fdc_size ! start after any existing fdc data
fdc_spt_coeffc_length = spt_size
fdc_spt_coeffs_start = fdc_spt_coeffc_start + fdc_spt_coeffc_length
fdc_spt_coeffs_length = spt_size

! calculate sizes and positions SKEB data in header array
fdc_skeb2_dpsidtc_start = fdc_spt_coeffs_start + fdc_spt_coeffs_length
fdc_skeb2_dpsidtc_length = skeb2_size
fdc_skeb2_dpsidts_start = fdc_skeb2_dpsidtc_start + fdc_skeb2_dpsidtc_length
fdc_skeb2_dpsidts_length = skeb2_size

! calculate sizes and positions of RP2 data in header array
fdc_rp2_coef_start = fdc_skeb2_dpsidts_start + fdc_skeb2_dpsidts_length
fdc_rp2_coef_length = rp_size

END SUBROUTINE calculate_fdc_positions

SUBROUTINE assign_a_flddepc_input_values()

IMPLICIT NONE

! Stochastic physics data stored in the dump header may have been created
! with different stochastic physics configurations / settings than requested
! in the current model run.  To avoid this causing problems in the
! stochastic physics schemes, the relevant data is first read in from the
! dump header to the array a_flddepc_in.  A new array, a_flddepc, is created
! and sized according to which stochastic physics schemes are specified at
! runtime.  Any relevant data from a_flddepc_in is then copied across to
! a_flddepc.  The new array a_flddepc is then used to store the stochastic
! physics data in the stochastic physics subroutines and is output to the
! header at dump time.

! In this routine, the array a_flddepc is initialised according to
! the stochastic physics schemes in use and the data available in
! a_flddepc_in.  If the size of a_flddepc is different from a_flddepc_in
! (such as may occur if different stochastic physics schemes or settings
! are used in the creation of the dump and the current forecast) then
! the fixed length header entries are updated with the new values.

INTEGER :: spt_size, skeb2_size, rp_size, fdc_size
INTEGER :: i, icode

CHARACTER(LEN=*), PARAMETER :: routinename = 'ASSIGN_A_FLDDEPC_INPUT_VALUES'

icode=0

! For each stph scheme that is switched on, check whether appropriately
! sized data exists in a_flddepc_in.  Set the size of data needed for
! each scheme based on namelist inputs (spt and skeb2) or parameter
! values (rp2). The size is set to zero if the scheme is not switched on.
! Update stph_header_flag to indicate which schemes will output data to
! the dump.

IF ( l_spt ) THEN
  IF ( MOD(a_inthd(ih_stochastic_flag)/stph_spt_data_present,2) == 1 ) THEN
    IF ( a_inthd(ih_stph_n2) == stph_n2 ) THEN
      stph_spt_data_check = stph_spt_data_present
    END IF
  END IF
  spt_size = (stph_n2 + 1) * stph_n2
  stph_header_flag = stph_header_flag + stph_spt_data_present
ELSE
  spt_size = 0
END IF

IF ( l_skeb2 ) THEN
  IF ( MOD(a_inthd(ih_stochastic_flag)/stph_skeb2_data_present,2) == 1 ) THEN
    IF ( a_inthd(ih_stph_n1) == stph_n1 .AND.                                 &
           a_inthd(ih_stph_n2) == stph_n2) THEN
      stph_skeb2_data_check = stph_skeb2_data_present
    END IF
  END IF
  skeb2_size = (stph_n2 + 1) * (stph_n2 - stph_n1+1)
  stph_header_flag = stph_header_flag + stph_skeb2_data_present
ELSE
  skeb2_size = 0
END IF

IF ( l_rp2 ) THEN
  IF ( MOD(a_inthd(ih_stochastic_flag)/stph_rp2_data_present,2) == 1 ) THEN
    stph_rp2_data_check = stph_rp2_data_present
  END IF
  rp_size = rp_max
  stph_header_flag = stph_header_flag + stph_rp2_data_present
ELSE
  rp_size = 0
END IF

! Calulate the position in the array for each scheme's data
fdc_size = 1
CALL calculate_fdc_positions( spt_size, skeb2_size, rp_size, fdc_size )

! total size to allocate in field of constants
fdc_size = fdc_size +                 &
           fdc_spt_coeffc_length +    &
           fdc_spt_coeffs_length +    &
           fdc_skeb2_dpsidtc_length + &
           fdc_skeb2_dpsidts_length + &
           fdc_rp2_coef_length

a_len1_flddepc = fdc_size
a_len2_flddepc = 1

IF (ALLOCATED(a_flddepc)) DEALLOCATE(a_flddepc)

icode=0
ALLOCATE( a_flddepc(a_len1_flddepc*a_len2_flddepc+1), STAT=icode,            &
          ERRMSG=umMessage )
IF (icode /= 0) THEN
  CALL ereport( modulename//routinename, icode                               &
                , 'Failure in allocating array [a_flddepc]'//      newline// &
                TRIM(umMessage) )
END IF

! Assign values to the array a_flddepc by checking each stph scheme
! in turn.  If the relevant data is present in the dump file then
! copy it from a_flddepc_in; otherwise set to rmdi.

a_flddepc(:) = rmdi

IF (l_spt) THEN
  IF ( stph_spt_data_check == stph_spt_data_present ) THEN
    WRITE(umMessage,'(A)') 'SPT data read in from dump header'
    CALL umPrint(umMessage,src='dump_headers_mod')
    DO  i = fdc_spt_coeffc_start, fdc_skeb2_dpsidtc_start - 1
      a_flddepc(i) = a_flddepc_in(i)
    END DO
  END IF
END IF

IF (l_skeb2) THEN
  IF ( stph_skeb2_data_check == stph_skeb2_data_present ) THEN
    WRITE(umMessage,'(A)') 'SKEB2 data read in from dump header'
    CALL umPrint(umMessage,src='dump_headers_mod')
    DO  i = fdc_skeb2_dpsidtc_start, fdc_rp2_coef_start - 1
      a_flddepc(i) = a_flddepc_in(i)
    END DO
  END IF
END IF

IF (l_rp2) THEN
  IF ( stph_rp2_data_check == stph_rp2_data_present ) THEN
    WRITE(umMessage,'(A)') 'RP2 data read in from dump header'
    CALL umPrint(umMessage,src='dump_headers_mod')
    DO  i = fdc_rp2_coef_start+1, fdc_rp2_coef_start + rp_size
      a_flddepc(i) = a_flddepc_in(i)
    END DO
  END IF
END IF

DEALLOCATE(a_flddepc_in)

! Check whether the size of a_flddepc is different from the data
! stored in the dump.  If it is, update the fixed length header
! entries accordingly

IF (a_fixhd(126) /= a_len1_flddepc) THEN

  ! Update the start position of dump elements that come after
  ! field dependent constants inline with new fdc_size

  ! compressed field 1 start
  IF (a_fixhd(140) > 0) THEN
    a_fixhd(140) = a_fixhd(140) - a_fixhd(126) + fdc_size
  END IF

  ! compressed field 2 start
  IF (a_fixhd(142) > 0) THEN
    a_fixhd(142) = a_fixhd(142) - a_fixhd(126) + fdc_size
  END IF

  ! compressed field 3 start
  IF (a_fixhd(144) > 0) THEN
    a_fixhd(144) = a_fixhd(144) - a_fixhd(126) + fdc_size
  END IF

  ! lookup start
  IF (a_fixhd(150) > 0) THEN
    a_fixhd(150) = a_fixhd(150) - a_fixhd(126) + fdc_size
  END IF

  ! data start
  IF (a_fixhd(160) > 0) THEN
    a_fixhd(160) = a_fixhd(160) - a_fixhd(126) + fdc_size
  END IF

  ! Finally store the size in the appropriate fixed header entry
  a_fixhd(126) = a_len1_flddepc
  a_fixhd(127) = a_len2_flddepc

END IF

END SUBROUTINE assign_a_flddepc_input_values

INTEGER FUNCTION calculate_fdc_start ()
! Calculate the start position of the additional parameters, based on start 
! positions and sizes of preceding header components

IMPLICIT NONE
INTEGER :: start_pos
  
  ! start of level dependant constants, as we know fixed, real, integer and 
  ! level-dependant constants must be present
  start_pos = a_fixhd(110) 
  IF (a_len1_levdepc > 0 .AND. a_len2_levdepc > 0) THEN
    start_pos = start_pos + a_len1_levdepc * a_len2_levdepc
  END IF 
  
  IF (a_len1_rowdepc > 0 .AND. a_len2_rowdepc > 0) THEN
    start_pos = start_pos + a_len1_rowdepc * a_len2_rowdepc
  END IF 

  IF (a_len1_coldepc > 0 .AND. a_len2_coldepc > 0) THEN
    start_pos = start_pos + a_len1_coldepc * a_len2_coldepc
  END IF 

  calculate_fdc_start = start_pos

END FUNCTION calculate_fdc_start

END MODULE dump_headers_mod
