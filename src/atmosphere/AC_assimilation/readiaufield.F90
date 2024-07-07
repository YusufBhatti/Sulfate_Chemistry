! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Read in and check a single IAU increment field

MODULE readiaufield_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='READIAUFIELD_MOD'

CONTAINS

SUBROUTINE ReadIAUField (                   &
                          IAU_unit,         & ! in
                          FixHd,            & ! in
                          Len1Lookup,       & ! in
                          Len2Lookup,       & ! in
                          Lookup,           & ! inout
                          IncNum,           & ! in
                          FieldNum,         & ! in
                          LocFldLen,        & ! in
                          FirstCallThisInc, & ! in
                          Field )             ! out

! Method:
!
!   1. Check field dimensions for compatibility with the model.
!   2. Read in field from IAU increment file.
!   3. If required, reset polar rows to their mean values.
!   4. If required, write basic field stats to standard output.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE IAU_mod, ONLY:   &
    IAU_NumFldCodes,  &
    IAU_FldCodes,     &
    IAU_FldDescs,     &
    L_IAU_IncDiags,   &
    L_IAU_ResetPoles

USE global_2d_sums_mod, ONLY: &
    global_2d_sums

USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParParams, ONLY: halo_type_no_halo, pnorth, psouth
USE Control_Max_Sizes
USE lookup_addresses
USE fieldstats_mod, ONLY: fieldstats
USE cppxref_mod,  ONLY: ppx_grid_type, ppx_atm_tall, ppx_halo_type
USE ppxlook_mod,  ONLY: exppxi, set_ppxi
USE umPrintMgr,   ONLY: umPrint, umMessage
USE nlsizes_namelist_mod, ONLY:                                    &
    global_land_field, global_row_length, global_rows, land_field, &
    len_fixhd, row_length, rows

USE model_domain_mod, ONLY: model_type, mt_global

USE readflds_mod, ONLY: readflds

IMPLICIT NONE


! Subroutine arguments:

INTEGER, INTENT(IN)    :: IAU_unit         ! File unit to read from
INTEGER, INTENT(IN)    :: FixHd(Len_FixHd) ! Fixed-length header
INTEGER, INTENT(IN)    :: Len1Lookup       ! First  dimension of lookup header
INTEGER, INTENT(IN)    :: Len2Lookup       ! Second dimension of lookup header

INTEGER, INTENT(INOUT) :: Lookup(Len1Lookup,Len2Lookup) ! Lookup header

INTEGER, INTENT(IN)    :: IncNum           ! Increment file number
INTEGER, INTENT(IN)    :: FieldNum         ! Field number in IAU increment file
INTEGER, INTENT(IN)    :: LocFldLen        ! Local length of field
LOGICAL, INTENT(IN)    :: FirstCallThisInc ! Is this the first call for this
                                           ! increment file?
REAL,    INTENT(OUT)   :: Field(LocFldLen) ! Local part of field

! Local constants:

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'READIAUFIELD'

! Local variables:

INTEGER :: i
INTEGER :: isec
INTEGER :: item
INTEGER :: Code
INTEGER :: Code_orig
INTEGER :: Level
INTEGER :: BufLen
INTEGER :: Inc_rows
INTEGER :: Inc_cols
INTEGER :: Model_rows
INTEGER :: Model_cols
INTEGER :: Model_size
INTEGER :: ICode
INTEGER :: IM_index
INTEGER :: grid_type
INTEGER :: halo_type_orig
INTEGER :: s_addr
INTEGER :: e_addr
INTEGER :: grid_type_orig

LOGICAL :: Code_changed

REAL :: polar_sum(1)
REAL :: polar_row(row_length, 1)
REAL :: polar_mean
REAL :: Global_max
REAL :: Global_min
REAL :: Global_mean
REAL :: Global_RMS

CHARACTER(LEN=800) :: readiaumessage
CHARACTER(LEN=7)   :: FldDesc
CHARACTER(LEN=4)   :: IncNumStr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!- End of header ---------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

Code   = Lookup(item_code, FieldNum)
Level  = Lookup(lblev,     FieldNum)
BufLen = Lookup(lblrec,    FieldNum)

!-------------------------------------------------------------------------------
! [1]: Check field dimensions for compatibility with the model.
!-------------------------------------------------------------------------------

Inc_rows = Lookup(lbrow, FieldNum)
Inc_cols = Lookup(lbnpt, FieldNum)

Model_rows = global_rows
Model_cols = global_row_length


IF (Code == 3) Model_rows = global_rows+1 ! v

Model_size = Model_rows * Model_cols

! Check field dimension compatibility with model
IF (Inc_rows == 0 .AND. Inc_cols == 0) THEN
  ! grid-type on Land/Sea points only

  IF (BufLen /= global_land_field) THEN

    ICode = 1
    readiaumessage(1:80)   = 'Field dimension mis-match.'
    readiaumessage(81:160) = ''
    WRITE (readiaumessage(161:240),*) '  Increment no.:     ', IncNum
    WRITE (readiaumessage(241:320),*) '  Field no.:         ', FieldNum
    WRITE (readiaumessage(321:400),*) '  STASH code:        ', Code
    WRITE (readiaumessage(401:480),*) '  Level:             ', Level
    WRITE (readiaumessage(481:560),*) '  Data land-points:  ', BufLen
    WRITE (readiaumessage(561:640),*) '  Model land-points: ', global_land_field
    WRITE (readiaumessage(641:800),*) ''

    CALL EReport (RoutineName, ICode, readiaumessage)

  END IF

ELSE

  IF (Inc_rows /= Model_rows .OR.  &
      Inc_cols /= Model_cols) THEN

    ICode = 1
    readiaumessage(1:80)   = 'Field dimension mis-match.'
    readiaumessage(81:160) = ''
    WRITE (readiaumessage(161:240),*) '  Increment no.:  ', IncNum
    WRITE (readiaumessage(241:320),*) '  Field no.:      ', FieldNum
    WRITE (readiaumessage(321:400),*) '  STASH code:     ', Code
    WRITE (readiaumessage(401:480),*) '  Level:          ', Level
    WRITE (readiaumessage(481:560),*) '  Increment rows: ', Inc_rows
    WRITE (readiaumessage(561:640),*) '  Increment cols: ', Inc_cols
    WRITE (readiaumessage(641:720),*) '  Model     rows: ', Model_rows
    WRITE (readiaumessage(721:800),*) '  Model     cols: ', Model_cols

    CALL EReport (RoutineName, ICode, readiaumessage)

  END IF
END IF

!-------------------------------------------------------------------------------
! [2]: Read in field from IAU increment file.
!-------------------------------------------------------------------------------

! readflds makes use of the PPXI data for the field being read in. However,
! this will only have been set up for qT if it has been requested as a
! STASH diagnostic. To get around this, we temporarily change its STASH code so
! that it is read in as if it were specific humidity.
IF (Code == 16207 .OR. Code == 18001) THEN
  Code_orig = Code
  Code      = 10
  Lookup(item_code, FieldNum) = Code
  Code_changed = .TRUE.
ELSE
  Code_changed = .FALSE.
END IF

! Temporarily change PPXI halo entry so that haloes are not read:
IM_index       = 1
isec           = Code / 1000
item           = MOD(Code, 1000)
halo_type_orig = exppxi(im_index, isec, item, ppx_halo_type)
CALL set_ppxi(isec, item, ppx_halo_type, halo_type_no_halo)

! When SMC increments are written as ancillaries (i.e. all theta points)
! we must temporarily adjust the grid-type in the PPXI array
grid_type_orig = exppxi(im_index, isec, item, ppx_grid_type)

IF (code == 9 .AND. buflen == Model_size) THEN

  CALL set_ppxi(isec, item, ppx_grid_type, ppx_atm_tall)

END IF

CALL readflds ( IAU_unit,   & ! in
                1,          & ! in
                FieldNum,   & ! in
                Lookup,     & ! in
                Field,      & ! out
                FixHd,      & ! in
                1,          & ! in
                ICode,      & ! out
                readiaumessage )    ! out

IF (ICode > 0) THEN
  WRITE(umMessage,*) 'ReadIAUField: Error reading IAU field no. ', FieldNum
  CALL umPrint(umMessage,src='readiaufield')
  CALL EReport (RoutineName, ICode, readiaumessage)
END IF

! Restore PPXI halo entry:
CALL set_ppxi(isec, item, ppx_halo_type, halo_type_orig)

! Restore STASH code:
IF (Code_changed) THEN
  Code = Code_orig
  Lookup(item_code, FieldNum) = Code
END IF

!-------------------------------------------------------------------------------
! [3]: If required, reset polar rows to their mean values.
!-------------------------------------------------------------------------------

grid_type = exppxi(1, isec, item, ppx_grid_type)
IF (L_IAU_ResetPoles          .AND. &
    model_type == mt_global .AND. &
    grid_type == ppx_atm_tall) THEN

  IF (at_extremity(PNorth)) THEN

    s_addr = 1 + row_length * (rows-1)
    e_addr = s_addr + row_length - 1

    polar_row(:,1) = Field(s_addr:e_addr)

    CALL global_2d_sums(polar_row, row_length, 1, 0, 0, 1, &
                        polar_sum, gc_proc_row_group)

    polar_mean = polar_sum(1) / REAL(global_row_length)

    Field(s_addr:e_addr) = polar_mean

  END IF ! (at_extremity(PNorth))

  IF (at_extremity(PSouth)) THEN

    s_addr = 1
    e_addr = s_addr + row_length - 1

    polar_row(:,1) = Field(s_addr:e_addr)

    CALL global_2d_sums(polar_row, row_length, 1, 0, 0, 1, &
                        polar_sum, gc_proc_row_group)

    polar_mean = polar_sum(1) / REAL(global_row_length)

    Field(s_addr:e_addr) = polar_mean

  END IF ! (at_extremity(PSouth))

END IF

!-------------------------------------------------------------------------------
! [4]: If required, write basic field stats to standard output.
!-------------------------------------------------------------------------------

IF (L_IAU_IncDiags) THEN

  IF (FirstCallThisInc) THEN

    ! Get increment number as a string:
    WRITE (IncNumStr,'(I4)') IncNum
    IncNumStr = ADJUSTL(IncNumStr)

    WRITE(umMessage,*) ''
    CALL umPrint(umMessage,src='readiaufield')
    WRITE(umMessage,'(A)') 'Summary of fields for IAU increment no.'// &
                           TRIM(IncNumStr) //':'
    CALL umPrint(umMessage,src='readiaufield')
    WRITE(umMessage,*) ''
    CALL umPrint(umMessage,src='readiaufield')
    WRITE(umMessage,*) '  Field   Level  Max          Min         ' &
                              //' Mean         RMS'
    CALL umPrint(umMessage,src='readiaufield')
    WRITE(umMessage,*) '  -----   -----  ---          ---         ' &
                              //' ----         ---'
    CALL umPrint(umMessage,src='readiaufield')
  END IF

  CALL FieldStats (                               &
                    LocFldLen,                    & ! in
                    Field,                        & ! in
                    grid_type,                    & ! in
                    halo_type_no_halo,            & ! in
                    Global_max,                   & ! out
                    Global_min,                   & ! out
                    Global_mean,                  & ! out
                    Global_RMS )                    ! out

  ! Get field description:
  DO i = 1, IAU_NumFldCodes
    IF (IAU_FldCodes(i) == Code) FldDesc = IAU_FldDescs(i)
  END DO

  WRITE(umMessage,'(3A,I4,A,4(A,ES12.5))')                                     &
    '   ', FldDesc, ' ', Level, ' ', ' ', Global_max,  ' ', Global_min, &
                                     ' ', Global_mean, ' ', Global_RMS
  CALL umPrint(umMessage,src='readiaufield')

END IF

! Restore PPXI grid_type entry:
CALL set_ppxi(isec, item, ppx_grid_type, grid_type_orig)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ReadIAUField
END MODULE readiaufield_mod
