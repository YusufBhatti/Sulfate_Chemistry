! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  write out fieldsfiles

MODULE crmstyle_write_ff_mask_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_WRITE_FF_MASK_MOD'

CONTAINS

SUBROUTINE crmstyle_write_ff_mask(it, ntimes)


USE crmstyle_cntl_mod, ONLY:                                             &
  num_x,num_y, mlevs

USE crmstyle_grid_info_mod, ONLY:                                        &
  mask_full_x, mask_full_y

USE crmwork_arrays_mod, ONLY:                                            &
  h_theta_sea

USE crmstyle_pp_data_mod, ONLY:                                         &
  iyear,imon,iday,ihour,imin,isec,isyear,ismon,isday,ishour,ismin,issec, &
  itim, bzy,bzx,bdy,bdx, mask_bzy, mask_bzx, mask_bdy, mask_bdx,         &
  pseudo_lat,pseudo_lon

USE missing_data_mod, ONLY: rmdi, imdi


USE crmstyle_output_hdr_mod, ONLY:      &
  bcu_mask_hdr                             ! UM Headers: bcu_mask_file

USE io

USE file_manager, ONLY: assign_file_unit

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

! subroutines

USE get_env_var_mod, ONLY: get_env_var

USE crmstyle_write_bcu_mask_mod,    ONLY: crmstyle_write_bcu_mask
USE ereport_mod, ONLY: ereport, ereport_finalise

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE packing_codes_mod, ONLY: PC_No_Packing

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Write out fieldsfiles of buoyant cloudy updraughts mask.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.5.

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER,INTENT(IN) ::    &
  it                     & ! call number
 ,ntimes                   ! Total number of times


!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k, ij,ii,jj                  ! loop counters

INTEGER ::               &
  icall_type             & ! type of call
 ,icall_type2              ! type of call
INTEGER ::               &
  ErrorStatus            & ! Error code from operations on file
 ,MaxFldsout               ! max output fields


CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_WRITE_FF_MASK"

TYPE(UM_Header_type)    :: example_hdr   ! Header to copy

TYPE(PP_Field_type) ::       &
  ref_modlev_field(mlevs)      ! pp fields on model levels
TYPE(PP_Field_type) ::       &
  ref_single_field(1)          ! pp fields - single level

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

! Set up reference pp header for output fields
! All output fields on same p grid (i.e. no uv grids etc)
! take date from input pp fields

! Setup output field headers for model level fields correctly

! integer header
ref_modlev_field(:) % Hdr % validyear  = iyear    ! take from input fields
ref_modlev_field(:) % Hdr % validmonth = imon
ref_modlev_field(:) % Hdr % validdate  = iday
ref_modlev_field(:) % Hdr % validhour  = ihour
ref_modlev_field(:) % Hdr % validmin   = imin
ref_modlev_field(:) % Hdr % validsec   = isec
ref_modlev_field(:) % Hdr % datayear  = isyear
ref_modlev_field(:) % Hdr % datamonth = ismon
ref_modlev_field(:) % Hdr % datadate  = isday
ref_modlev_field(:) % Hdr % datahour  = ishour
ref_modlev_field(:) % Hdr % datamin   = ismin
ref_modlev_field(:) % Hdr % datasec   = issec
!ref_modlev_field(:) % Hdr % lbtim     = 11    ! forecast & validity time
ref_modlev_field(:) % Hdr % lbtim     = itim    ! forecast & validity time
ref_modlev_field(:) % Hdr % fcrange   = 0      ! hours from input fields

ref_modlev_field(:) % Hdr % LBLRec     = mask_full_y*mask_full_x

ref_modlev_field(:) % Hdr % lbcode     = 101     ! rotated pole
ref_modlev_field(:) % Hdr % lbhem      = 3       ! limited area
ref_modlev_field(:) % Hdr % NumRows    = mask_full_y
ref_modlev_field(:) % Hdr % NumCols    = mask_full_x
ref_modlev_field(:) % Hdr % lbext      = 0       ! no extra data
ref_modlev_field(:) % Hdr % lbpack     = PC_No_Packing       ! unpacked output
ref_modlev_field(:) % Hdr % lbrel      = 3
ref_modlev_field(:) % Hdr % ppcode     = 0       ! none
ref_modlev_field(:) % Hdr % lbcfc      = 0       ! none
ref_modlev_field(:) % Hdr % lbproc     = 0

ref_modlev_field(:) % Hdr % lbvc       = 1       ! height
ref_modlev_field(:) % Hdr % lbrvc      = 0       ! none
ref_modlev_field(:) % Hdr % lbexp      = 0       ! none
ref_modlev_field(:) % Hdr % MO8proj    = 0       ! none
ref_modlev_field(:) % Hdr % MO8type    = 0       ! none
DO k=1,mlevs
  ref_modlev_field(k) % Hdr % MO8Level   = k       ! level number
END DO
DO k=1,4
  ref_modlev_field(:) % Hdr % LBRsvd(k) = 0    ! No values
END DO
ref_modlev_field(:) % Hdr % stcode   = 0       ! need to set for each fiels
ref_modlev_field(:) % Hdr % LBSrce   = 1111    ! reals
ref_modlev_field(:) % Hdr % LBUser1  = 1       ! reals
ref_modlev_field(:) % Hdr % LBUser3  = 0       !
ref_modlev_field(:) % Hdr % LBUser5  = 0       !
ref_modlev_field(:) % Hdr % LBUser6  = 0       !
ref_modlev_field(:) % Hdr % LBUser7  = 1       ! atmosphere model

! real part
DO k=1,mlevs
  ref_modlev_field(k) % Hdr %RLevel   = h_theta_sea(k)  ! height of level
END DO
ref_modlev_field(:) % Hdr % pseudolat  = pseudo_lat   ! from grid
ref_modlev_field(:) % Hdr % pseudolon  = pseudo_lon
ref_modlev_field(:) % Hdr % bgor       = 0.0
ref_modlev_field(:) % Hdr % ZerothLat  = mask_bzy   ! Zeroth, so deduct dy
ref_modlev_field(:) % Hdr % LatInt     = mask_bdy
ref_modlev_field(:) % Hdr % ZerothLon  = mask_bzx
ref_modlev_field(:) % Hdr % LonInt     = mask_bdx
ref_modlev_field(:) % Hdr % bmdi       = rmdi
ref_modlev_field(:) % Hdr % bmks       = 1.0

! allocate space for data
DO k=1,mlevs
  ALLOCATE ( ref_modlev_field(k) % RData(mask_full_y*mask_full_x ,1 ) )
END DO

!------------------------------------------------------------------------------
! Set up reference pp header for single fields
! take date from input pp fields 
!------------------------------------------------------------------------------
! integer header
ref_single_field(1) % Hdr % validyear  = iyear    ! take from input fields
ref_single_field(1) % Hdr % validmonth = imon
ref_single_field(1) % Hdr % validdate  = iday
ref_single_field(1) % Hdr % validhour  = ihour
ref_single_field(1) % Hdr % validmin   = imin
ref_single_field(1) % Hdr % validsec   = isec
ref_single_field(1) % Hdr % datayear  = isyear
ref_single_field(1) % Hdr % datamonth = ismon
ref_single_field(1) % Hdr % datadate  = isday
ref_single_field(1) % Hdr % datahour  = ishour
ref_single_field(1) % Hdr % datamin   = ismin
ref_single_field(1) % Hdr % datasec   = issec
! Now takes value from input so correct calendar info
ref_single_field(1) % Hdr % lbtim     = itim   ! forecast & validity time
ref_single_field(1) % Hdr % fcrange   = 0      ! hours from input fields

ref_single_field(1) % Hdr % LBLRec     = mask_full_y*mask_full_x

ref_single_field(1) % Hdr % lbcode     = 101     ! rotated pole
ref_single_field(1) % Hdr % lbhem      = 3       ! limited area
ref_single_field(1) % Hdr % NumRows    = mask_full_y
ref_single_field(1) % Hdr % NumCols    = mask_full_x
ref_single_field(1) % Hdr % lbext      = 0       ! no extra data
ref_single_field(1) % Hdr % lbpack     = PC_No_Packing       ! unpacked output
ref_single_field(1) % Hdr % lbrel      = 3
ref_single_field(1) % Hdr % ppcode     = 0       ! none 
ref_single_field(1) % Hdr % lbcfc      = 0       ! none
ref_single_field(1) % Hdr % lbproc     = 0
ref_single_field(1) % Hdr % lbexp      = 0       ! none
ref_single_field(1) % Hdr % MO8proj    = 0       ! none
ref_single_field(1) % Hdr % MO8type    = 0       ! none

ref_single_field(1) % Hdr % lbvc       = 129     ! surface
ref_single_field(1) % Hdr % lbrvc      = 0       ! none
ref_single_field(1) % Hdr % LBRsvd(:)  = 0      ! No values
ref_single_field(1) % Hdr % stcode     = 0       ! need to set for each fiels
ref_single_field(1) % Hdr % LBSrce     = 1111    ! 
ref_single_field(1) % Hdr % LBUser1    = 1       ! 
ref_single_field(1) % Hdr % LBUser3    = 0       ! 
ref_single_field(1) % Hdr % LBUser5    = 0       ! 
ref_single_field(1) % Hdr % LBUser6    = 0       ! 
ref_single_field(1) % Hdr % LBUser7    = 1       ! atmosphere model

ref_single_field(1) % Hdr % pseudolat  = pseudo_lat  ! from input fields
ref_single_field(1) % Hdr % pseudolon  = pseudo_lon
ref_single_field(1) % Hdr % bgor       = 0.0  
ref_single_field(1) % Hdr % ZerothLat  = mask_bzy   ! Zeroth
ref_single_field(1) % Hdr % LatInt     = mask_bdy
ref_single_field(1) % Hdr % ZerothLon  = mask_bzx 
ref_single_field(1) % Hdr % LonInt     = mask_bdx
ref_single_field(1) % Hdr % bmdi       = rmdi
ref_single_field(1) % Hdr % bmks       = 1.0

! allocate space for data

ALLOCATE ( ref_single_field(1) % RData( mask_full_y*mask_full_x,1 ) )

!=============================================================================
! Setup values for Fields File header
!=============================================================================
ALLOCATE (example_hdr %FixHd(LenFixHd))  ! fix header
example_hdr %FixHd(:) = imdi             ! intialise

! Need to be set ?
! Lengths of integer and real headers seem to need to correspond to pp headers
example_hdr % LenIntC  = 46    ! for a FF ? ( 15 for ancillary)
example_hdr % LenRealC = 38    ! for a FF ? (6 ancillary)
example_hdr % Len1LevDepC = mlevs
example_hdr % Len2LevDepC = 8
example_hdr % Len1RowDepC = 0
example_hdr % Len2RowDepC = 0
example_hdr % Len1ColDepC = 0
example_hdr % Len2ColDepC = 0
example_hdr % Len1FldsOfC = 0
example_hdr % Len2FldsOfC = 0
example_hdr % LenExtraC   = 0
example_hdr % LenHistFile = 0
example_hdr % LenCompFldI1 = 0
example_hdr % LenCompFldI2 = 0
example_hdr % LenCompFldI3 = imdi
example_hdr % Len1Lookup = 64
example_hdr % Len2Lookup = 2    ! small lookup table as used as example

! Note need to allocate these arrays or problems in init_pp called by new_hdr
ALLOCATE (example_hdr % RowDepC(1))  !
example_hdr % RowDepC    = 0.0
ALLOCATE (example_hdr % ColDepC(1))  !
example_hdr % ColDepC    = 0.0

!Positions set by a call to init_pp from new_hdr
example_hdr % FixHd(1)   = 20   !  version
example_hdr % FixHd(2)   = 1    ! Atmosphere
example_hdr % FixHd(3)   = 1    ! hybrid vertical coordinate type
example_hdr % FixHd(4)   = 3    ! LAM no wrap

example_hdr % FixHd(5)   = 3     ! fields file (example come from ancillary)
example_hdr % FixHd(6)   = 0     ! no run identifier
example_hdr % FixHd(8)   = 1     ! gregorian calendar
example_hdr % FixHd(9)   = 3     ! Arakawa C grid
example_hdr % FixHd(10)   = 0     ! gregorian calendar
example_hdr % FixHd(12)   = 803     ! UM version 803

example_hdr % FixHd(21)   = iyear    ! First validity time in FF
example_hdr % FixHd(22)   = imon     !
example_hdr % FixHd(23)   = iday     !
example_hdr % FixHd(24)   = ihour    !
example_hdr % FixHd(25)   = imin     !
example_hdr % FixHd(26)   = isec     !
example_hdr % FixHd(27)   = 0        ! Should be day number


example_hdr % FixHd(28)   = iyear    ! last validity time in FF
example_hdr % FixHd(29)   = imon     !
example_hdr % FixHd(30)   = iday     !
example_hdr % FixHd(31)   = ihour    !
example_hdr % FixHd(32)   = imin     !
example_hdr % FixHd(33)   = isec     !
example_hdr % FixHd(34)   = 0        ! Should be day number

! Assume always the same
example_hdr % FixHd(100) = 256+1    ! start of integer constants
example_hdr % FixHd(101) = 46       ! for a FF  ( 15 for ancillary)

example_hdr % FixHd(112) = imdi     !
example_hdr % FixHd(115) = imdi     ! No row-dependent constants
example_hdr % FixHd(116) = imdi     !
example_hdr % FixHd(117) = imdi     !
example_hdr % FixHd(120) = imdi     ! No column-dependent constants
example_hdr % FixHd(121) = imdi     !
example_hdr % FixHd(122) = imdi     !
example_hdr % FixHd(125:149)=  imdi     ! Not using

example_hdr % FixHd(151) = 64           ! 1st dimension lookup table
example_hdr % FixHd(152) = 2   ! 2nd dimension lookup table

ALLOCATE (example_hdr %IntC(example_hdr % LenIntC))
example_hdr %IntC(:) = imdi          ! initialise
example_hdr % IntC(6)    = mask_full_x       ! Number of points E-W
example_hdr % IntC(7)    = mask_full_y       ! Number of points N-S
example_hdr % IntC(8)    = mlevs
example_hdr % IntC(9)    = mlevs

example_hdr % IntC(21)    = imdi
ALLOCATE (example_hdr %RealC(example_hdr % LenRealC))
example_hdr %RealC(:) = rmdi          ! initialise

example_hdr % RealC(1)   = mask_bdx           ! E-W grid spacing in degrees
example_hdr % RealC(2)   = mask_bdy           ! N-S grid spacing in degrees
example_hdr % RealC(3)   = mask_bzy+mask_bdy  ! Latitude of first position
example_hdr % RealC(4)   = mask_bzx+mask_bdx  ! Longitude of first position
example_hdr % RealC(5)   = pseudo_lat         ! From original grid
example_hdr % RealC(6)   = pseudo_lon         ! From original grid

ALLOCATE (example_hdr %LevDepC(8*mlevs+1))
example_hdr %LevDepC(:) = rmdi          ! initialise
!=============================================================================

!-----------------------------------------------------------------------
! work out type of call

IF (it == 1) THEN
  icall_type = 1    ! first call need to open and set up header
ELSE IF (it == ntimes) THEN
  icall_type = 3    ! Final call need to close file and write header
ELSE
  icall_type = 2    ! No need to open or close file
END IF
!
!-----------------------------------------------------------------------
! BCu mask
!-----------------------------------------------------------------------
IF (icall_type == 1) THEN
  ! Only assign a unit number and open on first call
  CALL get_env_var("BCU_MASK", bcu_mask_hdr % FileName)

  CALL assign_file_unit(bcu_mask_hdr % FileName, &
                      bcu_mask_hdr % UnitNum, handler="portio")
END IF

maxfldsout = mlevs*2*ntimes  ! more than enough only 1 model level field
                             ! plus one single level field

CALL  crmstyle_write_bcu_mask(icall_type, maxfldsout, example_hdr,         &
                            ref_modlev_field, ref_single_field,            &
                            bcu_mask_hdr, errorstatus )


! Cleanup space - is this needed and will it help memory?

DEALLOCATE(example_hdr%LevDepC)
DEALLOCATE(example_hdr%RealC)
DEALLOCATE(example_hdr%IntC)
DEALLOCATE(example_hdr%ColDepC)
DEALLOCATE(example_hdr%RowDepC)
DEALLOCATE(example_hdr%FixHd)

DEALLOCATE(ref_single_field(1)%RData )   ! large in size as a full field

DO k=1,mlevs
  DEALLOCATE(ref_modlev_field(k)%RData )   ! large in size as a full field
END DO

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_write_ff_mask

END MODULE crmstyle_write_ff_mask_mod
