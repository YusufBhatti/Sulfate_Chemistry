! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!Purpose:      Calculates latitude and longitude on equatorial
!              latitude-longitude (eq) grid used in regional
!              models from input arrays of latitude and
!              longitude on standard grid. Both input and output
!              latitudes and longitudes are in degrees.

PROGRAM main

IMPLICIT NONE

INTEGER :: iflag,icode, j, io
REAL :: RLAT(2),RLON(2),OLAT(2),OLON(2)
REAL :: TLAT(2),tlon(2)
REAL :: zoom, plon, plat,xshift,yshift 
CHARACTER*80  cname,name
 
READ(5,*)zoom,plon,plat,xshift,yshift,icode
      
IF(icode == 0)name='../../data/data_coarse'
IF(icode == 1)name='../../data/data_fine'
IF(icode == 2)name='../../data/data_latlon'
 
OPEN(11,FILE=name,Form='UNFORMATTED')

DO J=1,100000
      
  READ(11,IOSTAT=io)TLON(1),TLAT(1),TLON(2),TLAT(2)
  IF (io < 0 ) EXIT 
  CALL LLTOEQ (TLAT,TLON,OLAT,OLON,PLAT,PLON,2)

  OLON(1)=OLON(1)*ZOOM+xshift
  OLON(2)=OLON(2)*ZOOM+xshift
  OLAT(1)=-OLAT(1)*ZOOM+yshift
  OLAT(2)=-OLAT(2)*ZOOM+yshift
  
  IF(OLON(1) > 2.*xshift -100. .AND. OLON(2) < 100.)THEN
     OLON(1)=AMOD(OLON(1),2.*xshift)-2.*xshift
  END IF
  
  IF(OLON(1) < 100. .AND. OLON(2) > 2.*xshift -100.)THEN
     OLON(2)=AMOD(OLON(2),2.*xshift)-2.*xshift
  END IF

  IF(OLON(1) < -10. .OR. OLAT(1) < -10. .OR.   &
     OLON(2) < -10. .OR. OLAT(2) < -10. .OR.   &
     OLON(1) > 2.*xshift+10. .OR. OLAT(1) > 2.*yshift+10. .OR. &
     OLON(2) > 2.*xshift+10. .OR. OLAT(2) > 2.*yshift+10.) THEN
  ELSE
     IF(IFLAG == 1)THEN
       RLON(1)=OLON(1)
       RLAT(1)=OLAT(1)
       RLON(2)=OLON(2)
       RLAT(2)=OLAT(2)
       IFLAG=1-IFLAG
     ELSE
       IF(RLON(2) == OLON(1).AND.RLAT(2) == OLAT(1))THEN
         WRITE(6,'(6(F7.1))')RLON(1),RLAT(1),RLON(2),RLAT(2), &
                             OLON(2),OLAT(2)
       ELSE
         WRITE(6,'(6(F7.1))')RLON(1),RLAT(1),RLON(2),RLAT(2), &
                             RLON(2),RLAT(2)
         WRITE(6,'(6(F7.1))')OLON(1),OLAT(1),OLON(2),OLAT(2), &
                             OLON(2),OLAT(2)
       END IF
       IFLAG=1-IFLAG
     END IF
   END IF
      
END DO
STOP
END
