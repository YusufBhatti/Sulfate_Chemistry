! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Parallel RCF: Transform from global to local co-ordinates:

MODULE Rcf_Global_To_Local_Mod

!  Subroutine Global_To_Local_Subdomain : global subdomain boundaries
!                                         to local ones
!  Subroutine global_to_local_RC: converts global row,column co-ords to
!                                 processor co-ordinates plus local
!                                 co-ordinates within the processor.
!   Function Rcf_Get_Fld_Type : Determines P, U or V type of field
!
! Description:
!   Takes a global definition of a subdomain region (in terms of
!   model gridpoints) and translates it into local numbers.
!   This effectively means local co-ordinates of the region of the
!   subdomain which intersects with this processor's area.
!
! Method:
!   Use the datastart variable in PARVARS to see if the requested
!   subdomain intersects with this processor's area, if it does
!   then use datastart to convert to local co-ordinate and do a bit
!   of logic using MAX and MIN to ensure the local co-ordinates
!   actually lie within the local area  Then make any corrections
!   necessary to account for a subdomain which crosses over the
!   0 longitude line. Finally, if L_include_halos is set to
!   .TRUE. - include any relevant halo regions.
!
! Derived from UM 4.5 code
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

USE field_types
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GLOBAL_TO_LOCAL_MOD'

CONTAINS

SUBROUTINE rcf_global_to_local_subdomain(                           &
                                   L_include_halosEW,               &
                                   L_include_halosNS,               &
                                   grid_code,procid,                &
                                   global_north_in,global_east_in,  &
                                   global_south_in,global_west_in,  &
                                   local_north,local_east,          &
                                   local_south,local_west)
USE UM_ParVars

USE UM_ParParams, ONLY: &
    halo_type_single

IMPLICIT NONE


! Subroutine arguments:
LOGICAL, INTENT(IN)   :: L_include_halosEW  ! If true, include East-West
                                            ! halos in local region
LOGICAL, INTENT(IN)   :: L_include_halosNS  ! if true incl. North-South
                                            ! halos in local region
INTEGER, INTENT(IN)   :: grid_code          ! STASH grid type
INTEGER, INTENT(IN)   :: procid             ! process result wanted for
INTEGER, INTENT(IN)   :: global_north_in    ! global northern boundary
INTEGER, INTENT(IN)   :: global_east_in     ! global eastern boundary
INTEGER, INTENT(IN)   :: global_south_in    ! global southern boundary
INTEGER, INTENT(IN)   :: global_west_in     ! global western boundary

INTEGER, INTENT(OUT)  :: local_north        ! local northern boundary
INTEGER, INTENT(OUT)  :: local_south        ! local sothern boundary
INTEGER, INTENT(OUT)  :: local_east         ! local eastern boundary
INTEGER, INTENT(OUT)  :: local_west         ! local westernboundary

INTEGER, PARAMETER    :: st_no_data = -3    ! Magic number

! Local variables
! Copies of the input arguments, that can be modified for
! wrap-around calculations
INTEGER   :: global_north, global_east, global_south, global_west
INTEGER   :: fld_type           ! is field on P or U or V grid?
INTEGER   :: row_len_nh         ! row length when halos are removed
INTEGER   :: nrows_nh           ! number of rows when halos are removed
INTEGER   :: first_global_pt_EW ! global point number of first and last
INTEGER   :: last_global_pt_EW  ! local points in local area
INTEGER   :: first_global_pt_NS ! in the East-West and
INTEGER   :: last_global_pt_NS  ! North-South directions

! Logicals indicating if this processor contains part of a
! subdomain
LOGICAL   :: NS_intersect, EW_intersect
LOGICAL   :: wrap      ! set to .TRUE. if the subdomain passes over the
                       ! the 0 degree longitude line
LOGICAL   :: fullfield ! if the field is NOT a subdomain
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_GLOBAL_TO_LOCAL_SUBDOMAIN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! ------------------------------------------------------------------

! Copy the global_in variables into local variables

global_north=global_north_in
global_east=global_east_in
global_south=global_south_in
global_west=global_west_in

! Find out if the data is on a mass or velocity grid

fld_type = Rcf_Get_Fld_Type(grid_code)

IF (fld_type  ==  fld_type_unknown) THEN
  WRITE(umMessage,*) 'GLOBAL_TO_LOCAL_SUBDOMAIN encountered ',     &
             'field with gridtype code ',grid_code
  CALL umPrint(umMessage,src='rcf_global_to_local_mod')
  WRITE(umMessage,*) 'Unable to process this field.'
  CALL umPrint(umMessage,src='rcf_global_to_local_mod')
  local_north=st_no_data
  local_south=st_no_data
  local_east=st_no_data
  local_west=st_no_data
  GO TO 9999
END IF

! Set up logical indicating if this is a full field, or just
! a subdomain

fullfield= ((( global_west  ==  1 ) .AND.            &
             ( global_south  ==  1 )) .AND.          &
           (((fld_type  ==  fld_type_p ) .AND.       &
             ( global_north  ==  glsizep(2) ) .AND.   &
             ( global_west   ==  glsizep(1) )) .OR.   &
            ((fld_type  ==  fld_type_u ) .AND.       &
             ( global_north  ==  glsizeu(2) ) .AND.  &
             ( global_west   ==  glsizeu(1))) .OR.   &
            ((fld_type  ==  fld_type_v ) .AND.       &
             ( global_north  ==  glsizev(2) ) .AND.  &
             ( global_west   ==  glsizev(1) ))))

! If this is a fullfield (ie. not a subdomain) the local addressing
! is easy:

IF (fullfield) THEN

  IF (L_include_halosNS) THEN
    local_north = g_lasize(2,fld_type_p, &
        halo_type_single,procid)
    local_south = 1
  ELSE
    local_north = g_lasize(2,fld_type_p, &
        halo_type_single,procid)-Offy
    local_south = 1+Offy
  END IF
  IF (L_include_halosEW) THEN
    local_west=1
    local_east=g_lasize(1,fld_type_p, &
        halo_type_single,procid)
  ELSE
    local_west=1+Offx
    local_east=g_lasize(1,fld_type_p, &
        halo_type_single,procid)-Offx
  END IF

ELSE ! a subdomain requires some careful analysis:

  IF (fld_type  ==  fld_type_p) THEN
    row_len_nh=g_blsizep(1,procid)
    nrows_nh=g_blsizep(2,procid)
  ELSE IF (fld_type  ==  fld_type_u) THEN
    row_len_nh=g_blsizeu(1,procid)
    nrows_nh=g_blsizeu(2,procid)
  ELSE
    row_len_nh=g_blsizev(1,procid)
    nrows_nh=g_blsizev(2,procid)
  END IF

  ! Set up variables giving the global point numbers of the
  ! start and end of this processor's subdomain

  first_global_pt_EW=g_datastart(1,procid)
  last_global_pt_EW=first_global_pt_EW+row_len_nh-1

  first_global_pt_NS=g_datastart(2,procid)
  last_global_pt_NS=first_global_pt_NS+nrows_nh-1

  ! If global_east is greater than the global row length, this
  ! indicates a wrap around - but not in the format this code
  ! wants - where it expects a wrap around to be indicated by
  ! the east column being less than the west column.

  IF (global_east  <   global_west) THEN
    wrap=.TRUE.
  ELSE IF (global_east  >   glsizep(1)) THEN
    wrap=.TRUE.
    global_east=global_east-glsizep(1)
  ELSE
    wrap=.FALSE.
  END IF

  EW_intersect =                                         &
         (( .NOT. wrap) .AND.                            &
          ((global_east  >=  first_global_pt_EW) .AND.   &
           (global_west  <=  last_global_pt_EW)))        &
         .OR.                                            &
         ((wrap) .AND.                                   &
          ((global_west  <=  last_global_pt_EW) .OR.     &
           (global_east  >=  first_global_pt_EW)))

  NS_intersect =                                         &
         ((global_south  <=  last_global_pt_NS) .AND.    &
          (global_north  >=  first_global_pt_NS))

  IF (NS_intersect) THEN

    IF ((global_south  >=  first_global_pt_NS) .AND.     &
        (global_south  <=  last_global_pt_NS)) THEN
      ! This processor contains the NS start of the subarea
      local_south=global_south-first_global_pt_NS+Offy+1
    ELSE
      ! This processor is to the North of the start of the subarea
      local_south=1+Offy
    END IF

    IF ((global_north  >=  first_global_pt_NS) .AND.     &
        (global_north  <=  last_global_pt_NS)) THEN
      ! This processor contains the NS end of the subarea
      local_north=global_north-first_global_pt_NS+Offy+1
    ELSE
      ! This processor is to the South of the subarea
      local_north=nrows_nh+Offy
    END IF

  ELSE

    local_north=st_no_data
    local_south=st_no_data

  END IF

  IF (EW_intersect) THEN

    IF ((global_west  >=  first_global_pt_EW) .AND.     &
        (global_west  <=  last_global_pt_EW)) THEN
      ! This processor contains the EW start of the subarea
      local_west=global_west-first_global_pt_EW+Offx+1
    ELSE
      ! This processor is to the right of the start of the subarea
      local_west=1+Offx
    END IF

    IF ((global_east  >=  first_global_pt_EW) .AND.     &
        (global_east  <=  last_global_pt_EW)) THEN
      ! This processor contains the EW end of the subarea
      local_east=global_east-first_global_pt_EW+Offx+1
    ELSE
      ! This processor is to the left of the end of the subarea
      local_east=Offx+row_len_nh
    END IF

  ELSE

    local_east=st_no_data
    local_west=st_no_data

  END IF

END IF ! is this a fullfield?

9999 CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Global_To_Local_Subdomain

!--------------------------------------------------------------------
! This routine is included for completeness
!--------------------------------------------------------------------

! Subroutine Interface:
SUBROUTINE Rcf_Global_To_Local_RC( grid_code,                      &
                                   global_column_in , global_row,  &
                                   processor_x , processor_y,      &
                                   local_column, local_row)

! Description:
! Takes a global co-ordinate, in model gridpoints, and returns
! the processor co-ordinate of the processor containing that
! point, and the local co-ordinates of the point on that processor.

USE UM_ParVars

IMPLICIT NONE


INTEGER, INTENT(IN)  :: grid_code          ! STASH grid type code
INTEGER, INTENT(IN)  :: global_column_in   ! global column number
INTEGER, INTENT(IN)  :: global_row         ! global row number
INTEGER, INTENT(OUT) :: processor_x        ! processor X (EW) co-ord.
                                           !       (0->nproc_x)
INTEGER, INTENT(OUT) :: processor_y        ! processor Y (NS) co-ord.
                                           ! (0->nproc_y)
INTEGER, INTENT(OUT) :: local_column       ! local column no. on proc.
INTEGER, INTENT(OUT) :: local_row          ! local row number on proc.

INTEGER, PARAMETER   :: st_no_data = -3


! Local variables

INTEGER   :: global_column ! modified version of global_column_in which
!                          ! takes account of domains wrapping over
!                          ! 0 degree longitude
INTEGER   :: fld_type      ! field stored on P grid or U grid?
INTEGER   :: row_len_nh,nrows_nh ! row_len and n_rows when halos removed
INTEGER   :: proc  ! loop counter for loop over processors

! global column and row numbers delimiting a processors area
INTEGER   :: start_col, end_col, start_row, end_row
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_GLOBAL_TO_LOCAL_RC'
! ------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Find out if the data is on a mass or velocity grid

fld_type=Rcf_Get_Fld_Type(grid_code)

IF (fld_type  ==  fld_type_unknown) THEN
  WRITE(umMessage,*) 'GLOBAL_TO_LOCAL_RC encountered ',       &
             'field with gridtype code ',grid_code
  CALL umPrint(umMessage,src='rcf_global_to_local_mod')
  WRITE(umMessage,*) 'Unable to process this field.'
  CALL umPrint(umMessage,src='rcf_global_to_local_mod')
  processor_x=st_no_data
  processor_y=st_no_data
  local_column=st_no_data
  local_row=st_no_data
  GO TO 9999
END IF

! If global_column_in is more than the global row length, perform
! a wrap around to ensure it falls within the global bounds

IF (global_column_in  >   glsizep(1)) THEN
  global_column=MOD(global_column_in+1,glsizep(1))-1
ELSE
  global_column=global_column_in
END IF

IF ((global_column  <   1) .OR.                        &
    (global_row  <   1) .OR.                           &
    (global_row  >   glsizep(2))) THEN

  WRITE(umMessage,*) 'GLOBAL_TO_LOCAL_RC encountered ',            &
             'impossible global row/column co-ordinates ', &
             'row: ',global_row,' column: ',global_column
  CALL umPrint(umMessage,src='rcf_global_to_local_mod')

  processor_x=st_no_data
  processor_y=st_no_data
  local_column=st_no_data
  local_row=st_no_data

END IF

! Make a first guess at the processor co-ordinates

processor_x=MIN(global_column/(glsizep(1)/gridsize(1)), nproc_x-1)
processor_y=MIN(global_row/(glsizep(2)/gridsize(2)), nproc_y-1)

proc=processor_x+processor_y*gridsize(1)

IF (fld_type  ==  fld_type_p) THEN
  row_len_nh=g_blsizep(1,proc)
  nrows_nh=g_blsizep(2,proc)
ELSE IF (fld_type  ==  fld_type_u) THEN
  row_len_nh=g_blsizeu(1,proc)
  nrows_nh=g_blsizeu(2,proc)
ELSE
  row_len_nh=g_blsizev(1,proc)
  nrows_nh=g_blsizev(2,proc)
END IF

start_col=g_datastart(1,proc)
end_col=start_col+row_len_nh-1
start_row=g_datastart(2,proc)
end_row=start_row+nrows_nh-1

! Now iterate around these processors until we hit the right one

DO WHILE  (((global_column  <   start_col) .OR.     &
            (global_column  >   end_col  ))         &
        .OR.                                        &
           ((global_row  <   start_row) .OR.        &
            (global_row  >   end_row)))


  IF (global_column  <   start_col) THEN
    processor_x=processor_x-1
  ELSE IF (global_column  >   end_col) THEN
    processor_x=processor_x+1
  END IF

  IF (global_row  <   start_row) THEN
    processor_y=processor_y-1
  ELSE IF (global_row  >   end_row) THEN
    processor_y=processor_y+1
  END IF

  proc=processor_x+processor_y*gridsize(1)

  IF (fld_type  ==  fld_type_p) THEN
    row_len_nh=g_blsizep(1,proc)
    nrows_nh=g_blsizep(2,proc)
  ELSE IF (fld_type  ==  fld_type_u) THEN
    row_len_nh=g_blsizeu(1,proc)
    nrows_nh=g_blsizeu(2,proc)
  ELSE
    row_len_nh=g_blsizev(1,proc)
    nrows_nh=g_blsizev(2,proc)
  END IF
  start_col=g_datastart(1,proc)
  end_col=start_col+row_len_nh-1
  start_row=g_datastart(2,proc)
  end_row=start_row+nrows_nh-1

END DO

! Now we have the processor co-ordinates, we can calculate the
! local co-ordinates.

local_column=Offx+global_column-start_col+1
local_row=Offy+global_row-start_row+1

9999 CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Global_To_Local_RC

!-------------------------------------------------------------------
!
! Rcf only version
!
!-------------------------------------------------------------------

! Function Interface
INTEGER FUNCTION rcf_get_fld_type (grid_type_code)

USE UM_ParVars
USE cppxref_mod, ONLY:                                            &
    ppx_atm_tall,       ppx_atm_tmerid,       ppx_atm_umerid,      &
    ppx_atm_uall,       ppx_atm_cuall,        ppx_atm_cvall,       &
    ppx_atm_river,      ppx_atm_ozone,        ppx_atm_compressed
IMPLICIT NONE


!
! Description:
! Takes a STASH grid type code, and returns which type of
! grid this is - mass or wind grid.

! Subroutine arguments:

INTEGER, INTENT(IN)   :: grid_type_code     ! IN : STASH grid type code

! Parameters

IF (((grid_type_code  >=  ppx_atm_tall)    .AND.    &
     (grid_type_code  <=  ppx_atm_tmerid)) .OR.     &
     (grid_type_code  ==  ppx_atm_ozone)   .OR.     &
     (grid_type_code  ==  ppx_atm_compressed)) THEN
  Rcf_Get_Fld_Type=fld_type_p
ELSE IF                                              &
   (((grid_type_code  >=  ppx_atm_uall)    .AND.    &
     (grid_type_code  <=  ppx_atm_umerid)) .OR.     &
     (grid_type_code  ==  ppx_atm_cuall))   THEN
  Rcf_Get_Fld_Type=fld_type_u
ELSE IF                                              &
   (grid_type_code  ==  ppx_atm_cvall)    THEN
  Rcf_Get_Fld_Type=fld_type_v
ELSE IF                                              &
   (grid_type_code  ==  ppx_atm_river) THEN
  Rcf_Get_Fld_Type=fld_type_r
ELSE
  Rcf_Get_Fld_Type=fld_type_unknown
END IF

RETURN

END FUNCTION Rcf_Get_Fld_Type

END MODULE Rcf_Global_To_Local_Mod



