! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routines: UNPACK --------------------------------------------------
!
! Purpose:  This subroutine uses the index arrays
!           INDEX_COMP(N_SEG), INDEX_EXP(N_SEG) and
!           INDEX_TO_ROWS(MAX_ROW,MAX_LEVEL) to unpack data
!           FROM T_COMP(N_COMP) into T_EXP(IMAX,JMAX,KMAX)
!
! Programming standard : UMDP PAPER 3
!
! -----------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE UNPACK(j1,j2,k1,k2,max_row,max_level,                  &
imax,jmax,kmax,index_comp,index_exp,n_seg,index_to_rows,          &
n_comp,t_comp,t_exp,real_mdi,cyclic)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!Argument variables
INTEGER ::                                                        &
 j1                                                               &
     ! (IN) FIRST ROW OF DATA TO BE UNPACKED
,j2                                                               &
     ! (IN) LAST ROW OF DATA TO BE UNPACKED
,k1                                                               &
     ! (IN) FIRST LEVEL OF DATA TO BE UNPACKED
,k2                                                               &
     ! (IN) LAST LEVEL OF DATA TO BE UNPACKED
,imax                                                             &
       ! (IN) NUMBER OF POINTS EAST-WEST
,jmax                                                             &
       ! (IN) NUMBER OF POINTS NORTH-SOUTH
,kmax                                                             &
       ! (IN) NUMBER OF POINTS IN VERTICAL
,n_seg                                                            &
       ! (IN) TOTAL NUMBER OF SEA SEGMENTS
,n_comp                                                           &
          ! (IN) DIMENSION OF COMPRESSED ARRAY
,max_row                                                          &
          ! (IN) UPPER LIMIT ON JMAX
,max_level                                                        &
            ! (IN) UPPER LIMIT ON KMAX
,index_comp(n_seg)                                                &
                    ! (IN) CONTAINS POSITIONS IN COMPRESSED ARRAY
!                         ! OF START OF EACH SEA SEGMENT
      ,index_exp(n_seg)                                                 &
                         ! (IN) CONTAINS POSITIONS IN 1-DIMENSIONAL
!                        ! EXPANDED ARRAY
!                        ! T_EXP_1_D(IMAX*MAX_ROW*MAX_LEVEL)
!                        ! OF START OF EACH SEA SEGMENT
      ,index_to_rows(max_row,max_level)                                 &
                                         ! (IN) CONTAINS NUMBER OF FIRST
!                                        ! SEA SEGMENT IN EACH ROW
!                                        ! AT EACH LEVEL IF THERE IS A
!                                        ! SEA SEGMENT IN THE ROW
!                                        ! CONTAINS NUMBER OF NEXT
!                                        ! SEA SEGMENT OTHERWISE
      ,i_data  ! NUMBER OF DISTINCT DATA POINTS EAST-WEST

REAL ::                                                           &
 real_mdi                                                         &
           ! (IN) REAL MISSING DATA INDICATOR
,t_comp(n_comp)                                                   &
                 ! (IN) COMPRESSED ARRAY
,t_exp(imax,jmax,kmax)  ! (OUT) 3-DIMENSIONAL EXPANDED ARRAY

LOGICAL ::                                                        &
 cyclic  ! (IN) INDICATES WHETHER T_EXP(IMAX,JMAX,KMAX) INCLUDES
!              ! DATA FOR CYCLIC WRAP-AROUND POINTS

!Local variables
INTEGER ::                                                        &
 j60                                                              &
              ! LOCAL LOOP INDEX FOR 60
,jspan                                                            &
              ! LOCAL NUMBER OF ITERATIONS OF ROWS
,kspan                                                            &
              ! LOCAL NUMBER OF ITERATIONS OF LEVELS
,i                                                                &
              ! LOCAL LOOP INDEX FOR COLUMNS
,j                                                                &
              ! LOCAL LOOP INDEX FOR ROWS
,k                                                                &
              ! LOCAL LOOP INDEX FOR LEVELS
,jn,kn                                                            &
              ! LOCAL VARIABLES
,num_seg                                                          &
              ! NUMBER OF SEA SEGMENTS IN PRESENT ROW
,len_seg                                                          &
              ! LENGTH OF PRESENT SEA SEGMENT
,COUNT                                                            &
              ! LOCAL COUNTER FOR POINTS IN A SEA SEGMENT
,seg                                                              &
              ! LOCAL LOOP INDEX FOR SEA SEGMENTS IN PRESENT ROW
,seg_pos                                                          &
              ! LOCAL COUNTER FOR SEA SEGMENTS
,x_pos                                                            &
              ! LOCAL COUNTER FOR POINTS IN A ROW
,ipoint_exp                                                       &
              ! LOCAL POINTER TO EXPANDED ARRAY
,ipoint_comp  ! LOCAL POINTER TO COMPRESSED ARRAY

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UNPACK'

!--------------------------------------------------------------------
!1. Fill the expanded array with real missing data indicators

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO k=1,kmax
  DO j=1,jmax
    DO i=1,imax
      t_exp(i,j,k)=real_mdi
    END DO
  END DO
END DO
!---------------------------------------------------------------------
!2. Set the wrap-around parameters

IF (cyclic) THEN
  i_data=imax-2
ELSE
  i_data=imax
END IF
!---------------------------------------------------------------------
!3. Unpack levels and rows
!3.1 Set the number of iterations over rows and levels
jspan = j2 - j1 + 1
kspan = k2 - k1 + 1

DO j60=0,jspan*kspan-1

  !3.3 Set up 2-d loop indices

  k = j60/jspan
  j = j60 - k*jspan + j1
  k = k + k1

  !3.4 Define the next row and level

  IF (j == max_row) THEN
    jn=1
    kn=k+1
  ELSE
    jn=j+1
    kn=k
  END IF

  !3.5 Calculate the number of segments in the present row

  IF (kn >  max_level) THEN
    num_seg=n_seg-index_to_rows(j,k)+1
  ELSE
    num_seg=index_to_rows(jn,kn)-index_to_rows(j,k)
  END IF

  DO seg=1,num_seg
    seg_pos=index_to_rows(j,k)+seg-1

    !3.6 Calculate the length of the present sea segment

    IF (seg_pos <  n_seg) THEN
      len_seg=index_comp(seg_pos+1)-index_comp(seg_pos)
    ELSE
      len_seg=n_comp-index_comp(seg_pos)+1
    END IF

    !3.7 Calculate t_exp(i,j,k) for each point in the segment

    DO COUNT=1,len_seg
      ipoint_exp=index_exp(seg_pos)+COUNT-1
      ipoint_comp=index_comp(seg_pos)+COUNT-1
      x_pos=ipoint_exp-(k-1)*i_data*max_row-(j-1)*i_data
      t_exp(x_pos,j-j1+1,k-k1+1)=t_comp(ipoint_comp)
    END DO
  END DO
END DO
!---------------------------------------------------------------------
!4. Calculate t_exp(i,j,k) for wrap-around points if necessary

IF (cyclic) THEN
  DO k=1,kmax
    DO j=1,jmax
      t_exp(imax-1,j,k)=t_exp(1,j,k)
      t_exp(imax,j,k)=t_exp(2,j,k)
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE UNPACK
