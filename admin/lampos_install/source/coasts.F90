! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

PROGRAM MAIN

IMPLICIT NONE

INTEGER :: j,k, io
REAL :: RLON1,RLAT1,RLON2,RLAT2
REAL :: OLON1,OLAT1,OLON2,OLAT2,OLON3,OLAT3,OLON4,OLAT4

OPEN(10,FILE='../data/data3')
OPEN(12,FILE='../data/data_coarse',FORM='UNFORMATTED')
OPEN(11,FILE='../data/data_fine',FORM='UNFORMATTED')
OPEN(13,FILE='../data/data_latlon',FORM='UNFORMATTED')

READ(10,'(2(f8.3))')RLON1,RLAT1
      
DO j=1,100000
  READ(10,'(2(f8.3))',IOSTAT=io)RLON2,RLAT2
  IF (io < 0 ) EXIT      
  IF(RLON2 == 0.)THEN
    READ(10,'(2(f8.3))',IOSTAT=io)RLON2,RLAT2
    IF (io < 0 ) EXIT
    RLAT1=RLAT2
    RLON1=RLON2
    IF(RLON1 == 0.)EXIT
  END IF

  OLON1=(RLON1-2.)*16. -180.
  OLAT1=(RLAT1-2.)*16. -90.
  OLON2=(RLON2-2.)*16. -180.
  OLAT2=(RLAT2-2.)*16. - 90.
  WRITE(11)OLON1,OLAT1,OLON2,OLAT2
  RLAT1=RLAT2
  RLON1=RLON2
      
END DO

REWIND 11

DO j=1,100000
  READ(11,IOSTAT=io)OLON1,OLAT1,OLON2,OLAT2
  IF (io < 0 ) EXIT
  READ(11,IOSTAT=io)OLON3,OLAT3,OLON4,OLAT4
  IF (io < 0 ) EXIT
  
  IF(OLON2 == OLON3 .AND. OLAT2 == OLAT3)THEN
    WRITE(12)OLON1,OLAT1,OLON4,OLAT4
  ELSE
    WRITE(12)OLON1,OLAT1,OLON2,OLAT2
    WRITE(12)OLON3,OLAT3,OLON4,OLAT4
  END IF  
END DO

! Lines of latitude longitude
DO k= -180,160,20
  DO j= 70,-65,-5
    WRITE(13)FLOAT(K)+.002,FLOAT(J)+.002,     &
             FLOAT(K)+.002,FLOAT(J-5)+.002
  END DO
END DO
      
DO j= 70,-70,-20
  DO k= -180,175,5
    WRITE(13)FLOAT(K)+.002,FLOAT(J)+.002,     &
             FLOAT(K+5)+.002,FLOAT(J)+.002 
  END DO
END DO

STOP
END
