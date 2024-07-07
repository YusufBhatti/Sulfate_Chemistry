! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Fills the LBC region of a LAM model with data from LBCs

SUBROUTINE set_lateral_boundaries(                                &
  row_length,rows,halo_i,halo_j,levels,fld_type,field,            &
  lenrim,lbc_size,lbc_start,lbc_halo_i,lbc_halo_j,lbc,            &
  rimwidth,n_rims_to_do,rimweights,at_extremity,                  &
  l_do_boundaries,l_do_halos)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParParams
IMPLICIT NONE

! Fills the LBC region of a LAM model FIELD with data from
! the LBC record.
! The N_RIMS_TO_DO variable allows only a certain depth of
! rim width to be done.
! The RIMWEIGHTS array contains a weighting factor that
! controls what proportion of FIELD is used and what
! proportion of LBC. In the halo area only LBC data is used.
!
! If RIMWIDTH is zero then the routine assumes there are no
! LBCs available and exits without updating any fields
!
! The logicals L_DO_BOUNDARIES and L_DO_HALOS allow each
! of the two components of the LBC region to be updated
! independently.
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
                    ! IN : size of FIELD halo in EW direction
, halo_j                                                          &
                    ! IN : size of FIELD halo in NS direction
, levels                                                          &
                    ! IN : number of vertical levels
, fld_type                                                        &
                    ! IN : type of input field (P,U or V)
, lenrim                                                          &
                    ! IN : size of one level of LBC data
, lbc_size(4)                                                     &
                    ! IN : size of each side (NESW) of LBC data
, lbc_start(4)                                                    &
                    ! IN : offset of each side of LBC data
, lbc_halo_i                                                      &
                    ! IN : size of LBC halo in EW direction
, lbc_halo_j                                                      &
                    ! IN : size of LBC halo in NS direction
, rimwidth                                                        &
                    ! IN : size (width) of the boundary area
, n_rims_to_do      ! IN : number of rims to do (counting from
                    !      the outside in)

LOGICAL ::                                                        &
  at_extremity(4)                                                 &
                    ! IN : Indicates if this processor is at
                    !      the edge (North,East,South,West)
                    !      of the processor grid
, l_do_boundaries                                                 &
                    ! IN : .TRUE. if the boundary region is
                    !      to be updated
, l_do_halos        ! IN : .TRUE. if the halo region is to
                    !      updated

REAL ::                                                           &
  rimweights(rimwidth)                                            &
                    ! IN : Weight to apply to each successive
                    !      rim
, lbc(lenrim,levels)                                              &
                    ! IN : LBC to apply to FIELD
, field(1-halo_i:row_length+halo_i,                               &
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
, lbc_row_len                                                     &
                    ! Length of a row of LBC data
, first_pt_of_lbc                                                 &
                    ! First point on row contained in LBC data
, first_row_of_lbc                                                &
                    ! First row contained in LBC data
, lbc_address                                                     &
                    ! address in the LBC array
, weights_I                                                       &
                    ! modified I to point to ew_weights array
, weights_J                                                       &
                    ! modified J to point to ns_weights array
, rim                                                             &
                    ! rim number (at corners)
, i,j,k                                                           &
                    ! loop counters (along row, row, level)
, rim_I             ! modified I to point to RIMWEIGHTS array

REAL ::                                                           &
  ns_weights(1-halo_i:row_length+halo_i,                          &
             1-halo_j:rimwidth)                                   &
, ew_weights(1-halo_i:rimwidth,                                   &
             1-halo_j:rows+halo_j)
                    ! Arrays which contain the weights
                    ! to apply to the LBCs at each point

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_LATERAL_BOUNDARIES'

!---------------------------------------------------------------
!
! The following diagram breaks up a LAM model area into a number
! of subcomponents (the letters will be referred to in the code):
! (Assumes the following model sizes for this example:
!  ROW_LENGTH=ROWS=12
!  RIMWIDTH=3
!  HALO_I=HALO_J=2)
!
!       North
!  aaaaaaaaaaaaaaaa
!  aaaaaaaaaaaaaaaa
!  bbcccddddddeeeff
!  bbcccddddddeeeff
!  bbcccddddddeeeff
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  gghhhiiiiiijjjkk
!  llmmmnnnnnnooopp
!  llmmmnnnnnnooopp
!  llmmmnnnnnnooopp
!  qqqqqqqqqqqqqqqq
!  qqqqqqqqqqqqqqqq
!       South

! Code

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!=================================================================
! 1. Set up the "*_weights" arrays

!-----------------------------------------------------------------
! 1.1 Do the North/South regions

IF (at_extremity(PNorth) .OR. at_extremity(PSouth)) THEN
  ! We will need the ns_weights array

  ! First we do the halo area (South : qlp, North : abf)

  ! The NS edge (South: q, North : a)
  DO j=1-halo_j,0
    DO i=1-halo_i,row_length+halo_i
      IF (l_do_halos) THEN
        ns_weights(i,j)=1.0
      ELSE ! don't update halos
        ns_weights(i,j)=0.0
      END IF ! IF (L_DO_HALOS
    END DO ! I
  END DO ! J

  ! The Western edge ( South: l, North : b)
  IF (at_extremity(PWest)) THEN
    DO j=1,rimwidth
      DO i=1-halo_i,0
        IF (l_do_halos) THEN
          ns_weights(i,j)=1.0
        ELSE ! don't update halos
          ns_weights(i,j)=0.0
        END IF ! IF (L_DO_HALOS)
      END DO ! I
    END DO ! J
  END IF ! IF (AT_EXTREMITY(PWest))

  ! The Eastern edge (South : p, North : f)
  IF (at_extremity(PEast)) THEN
    DO j=1,rimwidth
      DO i=row_length+1,row_length+halo_i
        IF (l_do_halos) THEN
          ns_weights(i,j)=1.0
        ELSE ! don't update halos
          ns_weights(i,j)=0.0
        END IF ! IF (L_DO_HALOS)
      END DO ! I
    END DO ! J
  END IF ! IF (AT_EXTREMITY(PEast))

  ! And now the boundary areas ( South : mno, North : cde )
  ! The Western corner (South : m, North : c)
  IF (at_extremity(PWest)) THEN
    DO j=1,rimwidth
      DO i=1,rimwidth
        IF (l_do_boundaries) THEN
          rim=MIN(i,j)
          IF (rim  <=  n_rims_to_do) THEN
            ns_weights(i,j)=rimweights(rim)
          ELSE ! This rim not required
            ns_weights(i,j)=0.0
          END IF !  IF (rim  <=  N_RIMS_TO_DO)
        ELSE ! don't update boundaries
          ns_weights(i,j)=0.0
        END IF ! IF (L_DO_BOUNDARIES)
      END DO ! I
    END DO ! J
  END IF ! IF (AT_EXTREMITY(PWest))

  ! The Eastern corner (South : o, North : e)
  IF (at_extremity(PEast)) THEN
    DO j=1,rimwidth
      DO i=row_length-rimwidth+1,row_length
        rim_I=row_length+1-i
        IF (l_do_boundaries) THEN
          rim=MIN(rim_I,j)
          IF (rim  <=  n_rims_to_do) THEN
            ns_weights(i,j)=rimweights(rim)
          ELSE ! This rim not required
            ns_weights(i,j)=0.0
          END IF !  IF (rim  <=  N_RIMS_TO_DO)
        ELSE ! don't update boundaries
          ns_weights(i,j)=0.0
        END IF ! IF (L_DO_BOUNDARIES)
      END DO ! I
    END DO ! J
  END IF ! IF (AT_EXTREMITY(PEast))

  ! The bit between the corners (South : n, North : d)
  IF (at_extremity(PEast)) THEN
    row_end_pt=row_length-rimwidth
  ELSE
    row_end_pt=row_length +halo_i
  END IF
  IF (at_extremity(Pwest)) THEN
    row_start_pt=rimwidth+1
  ELSE
    row_start_pt=1-halo_i
  END IF
  DO j=1,rimwidth
    DO i=row_start_pt,row_end_pt
      IF (l_do_boundaries) THEN
        IF (j  <=  n_rims_to_do) THEN
          ns_weights(i,j)=rimweights(j)
        ELSE ! This rim is not required
          ns_weights(i,j)=0.0
        END IF ! IF (J  <=  N_RIMS_TO_DO)
      ELSE ! don't update boundaries
        ns_weights(i,j)=0.0
      END IF ! IF (L_DO_BOUNDARIES)
    END DO ! I
  END DO ! J

END IF ! If we're at the North or South edge of the LAM

!-----------------------------------------------------------------
! 1.2 Do the East/West regions

IF (at_extremity(PWest) .OR. at_extremity(PEast)) THEN
  ! We will need the ew_weights array

  IF (at_extremity(PSouth)) THEN
    first_row=rimwidth+1
  ELSE ! not at the South
    first_row=1-halo_j
  END IF

  IF (at_extremity(PNorth)) THEN
    last_row=rows-rimwidth
  ELSE ! not at the North
    last_row=rows+halo_j
  END IF

  ! First the halo area (West : g , East k)

  DO j=first_row,last_row
    DO i=1-halo_i,0
      IF (l_do_halos) THEN
        ew_weights(i,j)=1.0
      ELSE ! don't update halos
        ew_weights(i,j)=0.0
      END IF ! IF (L_DO_HALOS)
    END DO ! I
  END DO ! J

  ! And now the boundary area (West : h, East : j)

  DO j=first_row,last_row
    DO i=1,rimwidth
      IF (l_do_boundaries) THEN
        IF (i  <=  n_rims_to_do) THEN
          ew_weights(i,j)=rimweights(i)
        ELSE ! this rim not required
          ew_weights(i,j)=0.0
        END IF ! IF (I  <=  N_RIMS_TO_DO)
      ELSE ! don't update boundaries
        ew_weights(i,j)=0.0
      END IF ! IF (L_DO_BOUNDARIES)
    END DO ! I
  END DO ! J

END IF ! IF at Western or Eastern edge

!=================================================================
! 2 Now apply the boundary conditions

!-----------------------------------------------------------------
! 2.1 Southern region

IF (at_extremity(PSouth)) THEN

  IF (l_do_halos) THEN
    first_row=1-halo_j
    row_start_pt=1-halo_i
    row_end_pt=row_length+halo_i
  ELSE
    first_row=1
    IF (at_extremity(PWest)) THEN
      row_start_pt = 1
    ELSE
      row_start_pt=1-halo_i
    END IF
    IF (at_extremity(PEast)) THEN
      row_end_pt = row_length
    ELSE
      row_end_pt = row_length+halo_i
    END IF
  END IF

  IF (l_do_boundaries) THEN
    last_row=rimwidth
  ELSE
    last_row=0
  END IF

  lbc_row_len=row_length + 2*lbc_halo_i

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( levels, first_row, last_row, row_start_pt, row_end_pt,   &
!$OMP         ns_weights, lbc_start, lbc_halo_j, lbc_row_len,          &
!$OMP         lbc_halo_i, field, lbc )                                 &
!$OMP PRIVATE( i, j, k, LBC_address )
  DO k=1,levels
    DO j=first_row,last_row
      DO i=row_start_pt,row_end_pt
        IF (ns_weights(i,j)  /=  0.0) THEN

          LBC_address=lbc_start(PSouth) +                         &
                      (j+lbc_halo_j-1)*lbc_row_len +              &
                      i+lbc_halo_i-1

          IF (ns_weights(i,j)  ==  1.0) THEN
            field(i,j,k)= lbc(LBC_address,k)
          ELSE
            field(i,j,k)=field(i,j,k)*(1.0-ns_weights(i,j)) +     &
                         lbc(LBC_address,k)*ns_weights(i,j)
          END IF
        END IF
      END DO ! I
    END DO ! J
  END DO ! K
!$OMP END PARALLEL DO

  IF ((.NOT. l_do_boundaries) .AND. l_do_halos) THEN
    first_row = 1
    last_row = rimwidth

    !      boundary area l
    IF (at_extremity(PWest)) THEN
      row_start_pt = 1-halo_i
      row_end_pt = 0
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( levels, first_row, last_row, row_start_pt, row_end_pt,   &
!$OMP         ns_weights, lbc_start, lbc_halo_j, lbc_row_len,          &
!$OMP         lbc_halo_i, field, lbc )                                 &
!$OMP PRIVATE( i, j, k, LBC_address )
      DO k=1,levels
        DO j=first_row,last_row
          DO i=row_start_pt,row_end_pt
            IF (ns_weights(i,j)  /=  0.0) THEN
              LBC_address=lbc_start(PSouth) +                     &
                      (j+lbc_halo_j-1)*lbc_row_len +              &
                      i+lbc_halo_i-1
              field(i,j,k)=lbc(LBC_address,k)
            END IF
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF   !  (AT_EXTREMITY(PWest))

    !      boundary area p
    IF (at_extremity(PEast)) THEN
      row_start_pt = row_length+1
      row_end_pt = row_length+halo_i

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( levels, first_row, last_row, row_start_pt, row_end_pt,   &
!$OMP         ns_weights, lbc_start, lbc_halo_j, lbc_row_len,          &
!$OMP         lbc_halo_i, field, lbc )                                 &
!$OMP PRIVATE( i, j, k, LBC_address )
      DO k=1,levels
        DO j=first_row,last_row
          DO i=row_start_pt,row_end_pt
            IF (ns_weights(i,j)  /=  0.0) THEN
              LBC_address=lbc_start(PSouth) +                     &
                      (j+lbc_halo_j-1)*lbc_row_len +              &
                      i+lbc_halo_i-1
              field(i,j,k)=lbc(LBC_address,k)
            END IF
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF    ! (AT_EXTREMITY(PEast))

  END IF      ! ((.not. L_DO_BOUNDARIES) .and. L_DO_HALOS)

END IF ! IF (AT_EXTREMITY(PSouth))

!-----------------------------------------------------------------
! 2.2 Northern region

IF (at_extremity(PNorth)) THEN

  IF (l_do_halos) THEN
    last_row=rows+halo_j
    row_start_pt=1-halo_i
    row_end_pt=row_length+halo_i
  ELSE
    last_row=rows
    IF (at_extremity(PWest)) THEN
      row_start_pt = 1
    ELSE
      row_start_pt=1-halo_i
    END IF
    IF (at_extremity(PEast)) THEN
      row_end_pt = row_length
    ELSE
      row_end_pt = row_length+halo_i
    END IF

  END IF

  IF (l_do_boundaries) THEN
    first_row=rows-rimwidth+1
  ELSE
    first_row=rows+1
  END IF

  lbc_row_len=row_length + 2*lbc_halo_i

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( levels, first_row, last_row, row_start_pt, row_end_pt,   &
!$OMP         ns_weights, lbc_start, rows, rimwidth, lbc_row_len,      &
!$OMP         lbc_halo_i, field, lbc )                                 &
!$OMP PRIVATE( i, j, k, weights_J, LBC_address )
  DO k=1,levels
    DO j=first_row,last_row
      weights_J=rows+1-j
      DO i=row_start_pt,row_end_pt
        IF (ns_weights(i,weights_J)  /=  0.0) THEN

          LBC_address=lbc_start(PNorth) +                         &
                      (j-(rows-rimwidth)-1)*lbc_row_len +         &
                      i+lbc_halo_i-1

          IF (ns_weights(i,weights_J)  ==  1.0) THEN
            field(i,j,k)=lbc(LBC_address,k)
          ELSE
            field(i,j,k)=                                         &
              field(i,j,k)*(1.0-ns_weights(i,weights_J)) +        &
              lbc(LBC_address,k)*ns_weights(i,weights_J)
          END IF
        END IF
      END DO ! I
    END DO ! J
  END DO ! K
!$OMP END PARALLEL DO

  IF ((.NOT. l_do_boundaries) .AND. l_do_halos) THEN
    first_row = rows-rimwidth+1
    last_row = rows

    !      boundary area b
    IF (at_extremity(PWest)) THEN
      row_start_pt = 1-halo_i
      row_end_pt = 0

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( levels, first_row, last_row, row_start_pt, row_end_pt,   &
!$OMP         ns_weights, lbc_start, rows, rimwidth, lbc_row_len,      &
!$OMP         lbc_halo_i, field, lbc )                                 &
!$OMP PRIVATE( i, j, k, weights_J, LBC_address )
      DO k=1,levels
        DO j=first_row,last_row
          weights_J=rows+1-j
          DO i=row_start_pt,row_end_pt
            IF (ns_weights(i,weights_J)  /=  0.0) THEN
              LBC_address=lbc_start(PNorth) +                     &
                      (j-(rows-rimwidth)-1)*lbc_row_len +         &
                      i+lbc_halo_i-1
              field(i,j,k)=lbc(LBC_address,k)
            END IF
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF   !  (AT_EXTREMITY(PWest))

    !      boundary area f
    IF (at_extremity(PEast)) THEN
      row_start_pt = row_length+1
      row_end_pt = row_length+halo_i

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( levels, first_row, last_row, row_start_pt, row_end_pt,   &
!$OMP         ns_weights, lbc_start, rows, rimwidth, lbc_row_len,      &
!$OMP         lbc_halo_i, field, lbc )                                 &
!$OMP PRIVATE( i, j, k, weights_J, LBC_address )
      DO k=1,levels
        DO j=first_row,last_row
          weights_J=rows+1-j
          DO i=row_start_pt,row_end_pt
            IF (ns_weights(i,weights_J)  /=  0.0) THEN
              LBC_address=lbc_start(PNorth) +                     &
                      (j-(rows-rimwidth)-1)*lbc_row_len +         &
                      i+lbc_halo_i-1
              field(i,j,k)=lbc(LBC_address,k)
            END IF
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF    ! (AT_EXTREMITY(PEast))

  END IF      ! ((.not. L_DO_BOUNDARIES) .and. L_DO_HALOS)

END IF !  IF (AT_EXTREMITY(PNorth))

!-----------------------------------------------------------------
! 2.3 Western region

IF (at_extremity(PWest)) THEN

  IF (l_do_halos) THEN
    row_start_pt=1-halo_i
  ELSE
    row_start_pt=1
  END IF

  IF (l_do_boundaries) THEN
    row_end_pt=n_rims_to_do
  ELSE
    row_end_pt=0
  END IF

  IF (at_extremity(PSouth)) THEN
    first_row=rimwidth+1
    first_row_of_LBC=rimwidth+1
  ELSE ! Not at the South
    first_row_of_LBC=1-lbc_halo_j
    first_row=1-halo_j
  END IF ! IF (AT_EXTREMITY(PSouth))

  IF (at_extremity(PNorth)) THEN
    last_row=rows-rimwidth
  ELSE ! Not at the North
    last_row=rows+halo_j
  END IF ! IF (AT_EXTREMITY(PNorth))

  lbc_row_len=lbc_halo_i+rimwidth

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( levels, first_row, last_row, row_start_pt, row_end_pt,   &
!$OMP         ew_weights, lbc_start, first_row_of_LBC, lbc_row_len,    &
!$OMP         lbc_halo_i, field, lbc )                                 &
!$OMP PRIVATE( i, j, k, LBC_address )
  DO k=1,levels
    DO j=first_row,last_row
      DO i=row_start_pt,row_end_pt
        IF (ew_weights(i,j)  /=  0.0) THEN

          LBC_address=lbc_start(PWest)+                           &
                      (j-first_row_of_LBC)*lbc_row_len +          &
                      i+lbc_halo_i-1

          IF (ew_weights(i,j)  ==  1.0) THEN
            field(i,j,k)=lbc(LBC_address,k)
          ELSE
            field(i,j,k)=field(i,j,k)*(1.0 - ew_weights(i,j)) +   &
                         lbc(LBC_address,k)*ew_weights(i,j)
          END IF
        END IF
      END DO ! I
    END DO ! J
  END DO ! 0K
!$OMP END PARALLEL DO

END IF ! IF (AT_EXTREMITY(PWest))

!-----------------------------------------------------------------
! 2.3 Eastern region

IF (at_extremity(PEast)) THEN

  IF (l_do_halos) THEN
    row_end_pt=row_length+halo_i
  ELSE
    row_end_pt=row_length
  END IF

  IF (l_do_boundaries) THEN
    row_start_pt=row_length-n_rims_to_do+1
  ELSE
    row_start_pt=row_length+1
  END IF

  IF (at_extremity(PSouth)) THEN
    first_row=rimwidth+1
    first_row_of_LBC=rimwidth+1
  ELSE ! Not at the South
    first_row_of_LBC=1-lbc_halo_j
    first_row=1-halo_j
  END IF ! IF (AT_EXTREMITY(PSouth))

  IF (at_extremity(PNorth)) THEN
    last_row=rows-rimwidth
  ELSE ! Not at the North
    last_row=rows+halo_j
  END IF ! IF (AT_EXTREMITY(PNorth))

  lbc_row_len=lbc_halo_i+rimwidth

  first_pt_of_LBC=row_length-rimwidth+1

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( levels, first_row, last_row, row_start_pt, row_end_pt,   &
!$OMP         ew_weights, lbc_start, first_row_of_LBC, lbc_row_len,    &
!$OMP         first_pt_of_LBC, field, lbc, rimwidth )                  &
!$OMP PRIVATE( i, j, k,  weights_I, LBC_address )
  DO k=1,levels
    DO j=first_row,last_row
      DO i=row_start_pt,row_end_pt
        weights_I=first_pt_of_LBC+rimwidth-i
        IF (ew_weights(weights_I,j)  /=  0.0) THEN

          LBC_address=lbc_start(PEast)+                           &
                      (j-first_row_of_LBC)*lbc_row_len +          &
                      i-first_pt_of_LBC

          IF (ew_weights(weights_I,j)  ==  1.0) THEN
            field(i,j,k)=lbc(LBC_address,k)
          ELSE
            field(i,j,k)=                                         &
              field(i,j,k)*(1.0 - ew_weights(weights_I,j)) +      &
              lbc(LBC_address,k)*ew_weights(weights_I,j)
          END IF
        END IF
      END DO ! I
    END DO ! J
  END DO ! K
!$OMP END PARALLEL DO

END IF ! IF (AT_EXTREMITY(PEast))

!=================================================================
! 3 Tidy up
!  Endgame:
!   If we're at the Eastern edge of a V-grid or P-grid field, the last point
!   on each row has not been updated. V-grid fields will be filled with zeros,
!   and P-grid fields with the value of their western neighbour.


! Zero last v for Endgame
!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( at_extremity, fld_type, levels, halo_j, rows, field,     &
!$OMP         row_length, halo_i )                                     &
!$OMP PRIVATE( j, k )
IF ( at_extremity(PEast) .AND.                                    &
     (fld_type  ==  fld_type_v)) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = 1-halo_j,rows+halo_j
      field(row_length+halo_i,j,k)=0.0
    END DO ! J
  END DO ! K
!$OMP END DO NOWAIT
END IF ! If at East and a V field (Endgame)

! Copy adjacent value for P field in Endgame
IF ( at_extremity(PEast) .AND.                                    &
     (fld_type  ==  fld_type_p)) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, levels
    DO j = 1-halo_j,rows+halo_j
      field(row_length+halo_i,j,k)=field(row_length+halo_i-1,j,k)
    END DO ! J
  END DO ! K
!$OMP END DO
END IF ! If at East and a P field (Endgame)
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_lateral_boundaries
