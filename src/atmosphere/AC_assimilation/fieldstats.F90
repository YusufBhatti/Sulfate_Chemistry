! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Calculate basic statistics for a horizontal field

MODULE fieldstats_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FIELDSTATS_MOD'

CONTAINS

SUBROUTINE FieldStats (              &
                        Field_len,   & ! in
                        Field,       & ! in
                        grid_type,   & ! in
                        halo_type,   & ! in
                        Global_max,  & ! out
                        Global_min,  & ! out
                        Global_mean, & ! out
                        Global_RMS )   ! out

! Description:
!
!   Calculate basic statistics for a horizontal field on u, v, theta or land
!   points.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE trignometric_mod, ONLY: &
    cos_theta_latitude,      &
    cos_v_latitude

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParParams, ONLY: halo_type_extended, halo_type_no_halo, halo_type_single
USE UM_ParCore, ONLY: nproc
USE Control_Max_Sizes
USE cppxref_mod, ONLY:       &
    ppx_atm_tall,             &
    ppx_atm_cuall,            &
    ppx_atm_cvall,            &
    ppx_atm_compressed

USE nlsizes_namelist_mod, ONLY: &
    land_field, n_rows, row_length, rows

USE missing_data_mod, ONLY: rmdi 

USE land_soil_dimensions_mod, ONLY: land_index

IMPLICIT NONE


! Subroutine arguments:

INTEGER, INTENT(IN)  :: Field_len
REAL,    INTENT(IN)  :: Field(Field_len) ! Horizontal field
INTEGER, INTENT(IN)  :: grid_type
INTEGER, INTENT(IN)  :: halo_type
REAL,    INTENT(OUT) :: Global_max
REAL,    INTENT(OUT) :: Global_min
REAL,    INTENT(OUT) :: Global_mean      ! Area-weighted mean
REAL,    INTENT(OUT) :: Global_RMS       ! Area-weighted RMS

! Local constants:

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'FIELDSTATS'

! Local variables:

INTEGER :: i, j, p
INTEGER :: ICode
INTEGER :: xm, ym
INTEGER :: xh, yh
INTEGER :: Expected_len

LOGICAL :: ThetaRows
LOGICAL :: VPoints
LOGICAL :: LandPoints

REAL :: val  
REAL :: Weight

REAL :: Global_sum
REAL :: Global_sumsq
REAL :: Global_sumwts

REAL :: ReshapedField(row_length,rows+1)
REAL :: WeightedField(row_length,rows+1)
REAL :: FieldWeights (row_length,rows+1)

REAL :: new_sum(3)

CHARACTER(LEN=320) :: fieldstatmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header ---------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! [1]: Do some checks based on grid and halo type.
!-------------------------------------------------------------------------------

ThetaRows  = .FALSE.
VPoints    = .FALSE.
LandPoints = .FALSE.

IF (grid_type == ppx_atm_tall .OR. &
    grid_type == ppx_atm_cuall) THEN
  ThetaRows = .TRUE.
  xm = row_length
  ym = rows
ELSE IF (grid_type == ppx_atm_cvall) THEN
  VPoints = .TRUE.
  xm = row_length
  ym = n_rows
ELSE IF (grid_type == ppx_atm_compressed) THEN
  LandPoints = .TRUE.
  Expected_len = field_len
ELSE
  ICode = 1
  fieldstatmessage = 'Grid type not catered for.'
  CALL EReport ( RoutineName, ICode, fieldstatmessage )
END IF

IF (.NOT. LandPoints) THEN

  IF (halo_type == halo_type_single) THEN
    xh = offx
    yh = offy
  ELSE IF (halo_type == halo_type_extended) THEN
    xh = halo_i
    yh = halo_j
  ELSE IF (halo_type == halo_type_no_halo) THEN
    xh = 0
    yh = 0
  ELSE
    ICode = 1
    fieldstatmessage = 'Invalid halo type.'
    CALL EReport ( RoutineName, ICode, fieldstatmessage )
  END IF

  Expected_len = (xm + 2*xh) * (ym + 2*yh)

END IF

IF (Field_len /= Expected_len) THEN
  ICode = 1
  fieldstatmessage(1:80)   = 'Supplied local field length incorrect.'
  fieldstatmessage(81:160) = ''
  WRITE (fieldstatmessage(161:240),*) '  supplied len:  ', Field_len
  WRITE (fieldstatmessage(241:320),*) '  expected len:  ', Expected_len
  CALL EReport ( RoutineName, ICode, fieldstatmessage )
END IF

!-------------------------------------------------------------------------------
! [2]: Get local statistics for fields on u, v, or theta points.
!-------------------------------------------------------------------------------

IF (.NOT. LandPoints) THEN
  global_max = -HUGE(1.0)
  global_min =  HUGE(1.0)
  new_sum(:) = 0.0

!$OMP PARALLEL DEFAULT(NONE)                                                 &
!$OMP PRIVATE(i,j)                                                           &
!$OMP SHARED(row_length,rows,halo_type,ym,xm,Field,xh,yh,cos_theta_latitude, &
!$OMP thetarows,cos_v_latitude,ReshapedField,WeightedField,FieldWeights)     &
!$OMP REDUCTION(MAX: global_max) REDUCTION(MIN: global_min)                  &
!$OMP REDUCTION(+: new_sum)

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      ReshapedField(i,j) = 0.0
      WeightedField(i,j) = 0.0
      FieldWeights (i,j) = 0.0
    END DO
  END DO
!$OMP END DO

  IF (halo_type == halo_type_no_halo) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, ym
      DO i = 1, xm
        ReshapedField(i,j) = Field(i + (j-1)*xm)
      END DO
    END DO
!$OMP END DO
  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, ym
      DO i = 1, xm
        ReshapedField(i,j) = Field(i+xh + (yh+j-1)*(xm+2*xh))
      END DO
    END DO
!$OMP END DO
  END IF

  IF (ThetaRows) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, ym
      DO i = 1, xm
        FieldWeights(i,j) = cos_theta_latitude(i,j)
      END DO
    END DO
!$OMP END DO
  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, ym
      DO i = 1, xm
        FieldWeights(i,j) = cos_v_latitude(i,j)
      END DO
    END DO
!$OMP END DO
  END IF

  ! Ignore undefined values in statistics
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, ym
    DO i = 1, xm
      IF (ReshapedField(i,j) == rmdi) FieldWeights(i,j) = 0.0
      global_min = MIN(global_min, ReshapedField(i,j) )
      global_max = MAX(global_max, ReshapedField(i,j) )
      WeightedField(i,j) = ReshapedField(i,j) &
                         * FieldWeights (i,j)
      new_sum(1) = new_sum(1) + WeightedField(i,j)
      new_sum(2) = new_sum(2) + ReshapedField(i,j) &
                              * WeightedField(i,j)
      new_sum(3) = new_sum(3) + FieldWeights (i,j)
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

END IF

!-------------------------------------------------------------------------------
! [3]: Get local statistics for fields on land points.
!-------------------------------------------------------------------------------

IF (LandPoints) THEN

  global_max = -HUGE(1.0)
  global_min =  HUGE(1.0)
  new_sum(:) = 0.0
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP PRIVATE(i,j,p,val,weight)                                       &
!$OMP SHARED(land_field,land_index,row_length,Field,cos_theta_latitude) &
!$OMP REDUCTION(MAX: global_max) REDUCTION(MIN: global_min)             &
!$OMP REDUCTION(+: new_sum)
  DO p = 1, land_field
    i = 1 + MOD(land_index(p)-1, row_length)
    j = 1 + (land_index(p)-1) / row_length
    val        = Field(p)
    Weight     = cos_theta_latitude(i,j)
    global_max = MAX(global_max, val)
    global_min = MIN(global_min, val)
    new_sum(1) = new_sum(1) + Weight * val  
    new_sum(2) = new_sum(2) + Weight * val**2
    new_sum(3) = new_sum(3) + Weight
  END DO
!$OMP END PARALLEL DO

END IF

!-------------------------------------------------------------------------------
! [4]: Get global statistics.
!-------------------------------------------------------------------------------

CALL gc_rmax(1, nproc, Icode, Global_max)
CALL gc_rmin(1, nproc, Icode, Global_min)
CALL gc_rsum(3, nproc, Icode, new_sum)

Global_mean =       new_sum(1) / new_sum(3)
Global_RMS  = SQRT( new_sum(2) / new_sum(3) )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE FieldStats
END MODULE fieldstats_mod
