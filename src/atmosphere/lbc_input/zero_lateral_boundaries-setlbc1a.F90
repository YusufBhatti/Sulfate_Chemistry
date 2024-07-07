! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Fills the LBC region of a LAM model with zeros

SUBROUTINE zero_lateral_boundaries(                               &
  row_length,rows,halo_i,halo_j,levels,fld_type,field,            &
  rimwidth,at_extremity,                                          &
  l_zero_boundaries,l_zero_halos)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParParams
IMPLICIT NONE

!  Fills the LBC region of a LAM model with zeros. The region can
!  be split into two areas:
!  1) Boundary - true edge of the model data
!  2) Halo - extended "halo" edge of the model data
!  Both of these areas can be zeroed independently
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: LBC Input
!
! Arguments:

INTEGER ::                                                        &
  row_length                                                      &
                    ! IN : number of points in row
, rows                                                            &
                    ! IN : number of rows
, halo_i                                                          &
                    ! IN : size of halo in EW direction
, halo_j                                                          &
                    ! IN : size of halo in NS direction
, levels                                                          &
                    ! IN : number of vertical levels
, fld_type                                                        &
                    ! IN : type of input field (P,U or V)
, rimwidth          ! IN : size (width) of the boundary area

LOGICAL ::                                                        &
  at_extremity(4)                                                 &
                    ! IN : Indicates if this processor is at
                    !      the edge (North,East,South,West)
                    !      of the processor grid
, l_zero_boundaries                                               &
                    ! IN : Do we want to zero the
                    !      boundary region
, l_zero_halos      ! IN : Do we want to zero the halo
                    !      region


REAL ::                                                           &
  field(1-halo_i:row_length+halo_i,                               &
                                     ! IN/OUT : field to
        1-halo_j:rows+halo_j,                                     &
                                     !          update
        levels)

! Local variables
INTEGER ::                                                        &
  row_start_pt                                                    &
                    ! first point along row to update
, row_end_pt                                                      &
                    ! last point along row to update
, first_row                                                       &
                    ! first row to update
, last_row                                                        &
                    ! last row to update
, i,j,k             ! loop counters (along row, row)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ZERO_LATERAL_BOUNDARIES'

! Code

!--------------------------------------------------------------------
! Northern region

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (at_extremity(PNorth)) THEN  ! This processor has a northern
                                ! LBC region
  IF (l_zero_halos) THEN        ! Zero out the halo area of the
                                ! LBC region
    row_start_pt=1-halo_i
    row_end_pt=row_length+halo_i
    first_row=rows+1
    last_row=rows+halo_j

    DO k=1,levels
      DO j=first_row,last_row
        DO i=row_start_pt,row_end_pt
          field(i,j,k)=0.0
        END DO
      END DO
    END DO

  END IF !  IF (L_ZERO_HALOS)

  IF (l_zero_boundaries) THEN   ! Zero out the boundary area
                                ! of the LBC region
    IF (at_extremity(PWest)) THEN
      row_start_pt = 1
    ELSE
      row_start_pt = 1-halo_i
    END IF

    IF (at_extremity(PEast)) THEN
      IF (fld_type  ==  fld_type_u) THEN
        row_end_pt=row_length-1
      ELSE
        row_end_pt=row_length
      END IF
    ELSE
      IF (fld_type  ==  fld_type_u) THEN
        row_end_pt=row_length+halo_i-1
      ELSE
        row_end_pt=row_length+halo_i
      END IF
    END IF

    first_row=rows-rimwidth+1
    last_row=rows

    DO k=1,levels
      DO j=first_row,last_row
        DO i=row_start_pt,row_end_pt
          field(i,j,k)=0.0
        END DO
      END DO
    END DO

  END IF ! IF (L_ZERO_BOUNDARIES)

END IF ! IF (AT_EXTREMITY(PNorth))

!--------------------------------------------------------------------
! Eastern region

IF (at_extremity(PEast)) THEN   ! This processor has a eastern
                                ! LBC region
  IF (l_zero_halos) THEN        ! Zero out the halo area of the
                                ! LBC region

    IF (fld_type  ==  fld_type_u) THEN
      row_start_pt=row_length
    ELSE
      row_start_pt=row_length+1
    END IF
    row_end_pt=row_length+halo_i

    IF (at_extremity(PSouth)) THEN
      first_row=1
    ELSE
      first_row=1-halo_j
    END IF
    IF (at_extremity(PNorth)) THEN
      last_row=rows
    ELSE
      last_row=rows+halo_j
    END IF

    DO k=1,levels
      DO j=first_row,last_row
        DO i=row_start_pt,row_end_pt
          field(i,j,k)=0.0
        END DO
      END DO
    END DO

  END IF !  IF (L_ZERO_HALOS)

  IF (l_zero_boundaries) THEN   ! Zero out the boundary area
                                ! of the LBC region
    IF (fld_type  ==  fld_type_u) THEN
      row_start_pt=row_length-rimwidth
      row_end_pt=row_length-1
    ELSE
      row_start_pt=row_length-rimwidth+1
      row_end_pt=row_length
    END IF
    IF (at_extremity(PSouth)) THEN
      first_row = 1
    ELSE
      first_row = 1-halo_j
    END IF

    IF (at_extremity(PNorth)) THEN
      last_row = rows
    ELSE
      last_row = rows+halo_j
    END IF


    DO k=1,levels
      DO j=first_row,last_row
        DO i=row_start_pt,row_end_pt
          field(i,j,k)=0.0
        END DO
      END DO
    END DO

  END IF ! IF (L_ZERO_BOUNDARIES)

END IF ! IF (AT_EXTREMITY(PEast))

!--------------------------------------------------------------------
! Southern region

IF (at_extremity(PSouth)) THEN   ! This processor has a eastern
                                ! LBC region
  IF (l_zero_halos) THEN        ! Zero out the halo area of the
                                ! LBC region
    row_start_pt=1-halo_i
    row_end_pt=row_length+halo_i
    first_row=1-halo_j
    last_row=0

    DO k=1,levels
      DO j=first_row,last_row
        DO i=row_start_pt,row_end_pt
          field(i,j,k)=0.0
        END DO
      END DO
    END DO

  END IF !  IF (L_ZERO_HALOS)

  IF (l_zero_boundaries) THEN   ! Zero out the boundary area
                                ! of the LBC region
    IF (at_extremity(PWest)) THEN
      row_start_pt = 1
    ELSE
      row_start_pt = 1-halo_i
    END IF

    IF (at_extremity(PEast)) THEN
      IF (fld_type  ==  fld_type_u) THEN
        row_end_pt=row_length-1
      ELSE
        row_end_pt=row_length
      END IF
    ELSE
      IF (fld_type  ==  fld_type_u) THEN
        row_end_pt=row_length+halo_i-1
      ELSE
        row_end_pt=row_length+halo_i
      END IF
    END IF

    first_row=1
    last_row=rimwidth

    DO k=1,levels
      DO j=first_row,last_row
        DO i=row_start_pt,row_end_pt
          field(i,j,k)=0.0
        END DO
      END DO
    END DO

  END IF ! IF (L_ZERO_BOUNDARIES)

END IF ! IF (AT_EXTREMITY(PNorth))

!--------------------------------------------------------------------
! Western region

IF (at_extremity(PWest)) THEN   ! This processor has a western
                                ! LBC region
  IF (l_zero_halos) THEN        ! Zero out the halo area of the
                                ! LBC region
    row_start_pt=1-halo_i
    row_end_pt=0

    IF (at_extremity(PSouth)) THEN
      first_row=1
    ELSE
      first_row=1-halo_j
    END IF
    IF (at_extremity(PNorth)) THEN
      last_row=rows
    ELSE
      last_row=rows+halo_j
    END IF

    DO k=1,levels
      DO j=first_row,last_row
        DO i=row_start_pt,row_end_pt
          field(i,j,k)=0.0
        END DO
      END DO
    END DO

  END IF !  IF (L_ZERO_HALOS)

  IF (l_zero_boundaries) THEN   ! Zero out the boundary area
                                ! of the LBC region
    row_start_pt=1
    row_end_pt=rimwidth
    IF (at_extremity(PSouth)) THEN
      first_row = 1
    ELSE
      first_row = 1-halo_j
    END IF

    IF (at_extremity(PNorth)) THEN
      last_row = rows
    ELSE
      last_row = rows+halo_j
    END IF


    DO k=1,levels
      DO j=first_row,last_row
        DO i=row_start_pt,row_end_pt
          field(i,j,k)=0.0
        END DO
      END DO
    END DO

  END IF ! IF (L_ZERO_BOUNDARIES)

END IF ! IF (AT_EXTREMITY(PEast))

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE zero_lateral_boundaries
