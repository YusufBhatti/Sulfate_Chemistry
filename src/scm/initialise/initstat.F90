! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine initstat
! Purpose:-           To calculate the initial variables required by
!                     statistical forcing routines used later
!                     and also prints out initial climate datasets
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!=====================================================================

SUBROUTINE initstat                                                           &
  ( row_length, rows, ntrop, andayy, dayno, q, t, lat, longitude              &
  , p_in, pa, pb, alfada, alfadb, tbara, tbarb, tsda, tsdb, tgrada, tgradb    &
  , dbara, dbarb, dgrada, dgradb, vnbara, vnbarb, vnsda, vnsdb, vpbara        &
  , vpbarb, wbara, wbarb, wsda, wsdb, atime, btime, p_theta_levels )

USE conversions_mod, ONLY: pi

USE nlsizes_namelist_mod, ONLY: model_levels

USE qsat_mod, ONLY: qsat

IMPLICIT NONE

INTEGER ::     &
  row_length   &! In x direction dimension
, rows         &! In y direction dimension
, ntrop         ! In Max number of levels in the troposphere

!---------------------------------------------------------------------
!     Arguments
!---------------------------------------------------------------------

REAL ::         &
  andayy         ! In No. of days in 1 year

INTEGER ::      &
  dayno          ! In  Day number relative to winter solstice

REAL ::                                 &
  alfada(row_length,rows)               &! Out Amplitude and mean of
, alfadb(row_length,rows)               &!     seasonal variation of tuning
                                         !     hybrid vertical coordinate
, atime,btime                           &! Out Constants for calculating
                                         !     annual cycle used in eqn 2.33
                                         !     in SCM doc.
, dbara(row_length,rows,model_levels)   &! Out Amplitude and mean of seasonal
, dbarb(row_length,rows,model_levels)   &!     variation of mean dew pt.
                                         !     depression (K)
, dgrada(row_length,rows,model_levels)  &! Out Amplitude and mean of seasonal
, dgradb(row_length,rows,model_levels)  &!     variation of gradient of
                                         !     dew pt. depression (K/km)
, lat0                                  &! Dummy for I/Os
, lat(row_length,rows)                  &! Out Latitude and longitude of
, longitude(row_length,rows)            &!     gridpoint
, longitude0                            &! Dummy for I/Os
, pa(row_length,rows, model_levels+1)   &! Out Amplitude and mean of seasonal
, pb(row_length,rows, model_levels+1)   &!     variation of pressure
, q(row_length,rows,model_levels)       &! InOut Specific humidity (kg/kg)
, p_theta_levels(row_length,rows,model_levels) &
                                         ! In pressure ((HPa or mb))
, t(row_length,rows,model_levels)       &! InOut Temps(K)
, tbara(row_length,rows,model_levels)   &! Out Amplitude and mean of seasonal
, tbarb(row_length,rows,model_levels)   &!     variation of mean temo. (K)
, tgrada(row_length,rows,model_levels)  &! Out Amplitude and mean of seasonal
, tgradb(row_length,rows,model_levels)  &!    variation of temp. gradient (K/km)
, tsda(row_length,rows,model_levels)    &! Out Amplitude and mean of seasonal
, tsdb(row_length,rows,model_levels)    &!     variation of SD of temp. (K)
, tstara(row_length,rows)               &! Out Amplitude and mean of seasonal
, tstarb(row_length,rows)               &!     variation of surface temp. (K)
, vnbara(row_length,rows,model_levels)  &! Out Amplitude and mean of seasonal
, vnbarb(row_length,rows,model_levels)  &!     variation of velocity VN (m/s)
, vnsda(row_length,rows,model_levels)   &! Out Amplitude and mean of seasonal
, vnsdb(row_length,rows,model_levels)   &!  variation of SD of velocity VN (m/s)
, vpbara(row_length,rows,model_levels)  &! Out Amplitude and mean of seasonal
, vpbarb(row_length,rows,model_levels)  &!     variation of velocity VP (m/s)
, wbara(row_length,rows,ntrop)          &! Out Amplitude and mean of seasonal
, wbarb(row_length,rows,ntrop)          &!     variation of vert. vel.( mb/s)
, wsda(row_length,rows,ntrop)           &! Out Amplitude and mean of seasonal
, wsdb(row_length,rows,ntrop)            !  variation of SD of vert. vel. (mb/s)

!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------

! Variables for printout of climate dataset

CHARACTER(LEN=29) ::   &
  cfmt             ! Format statement for each row of variables

CHARACTER(LEN=34) ::   &
  ctfmt            ! Format statement for title of each row

INTEGER ::              &
  element               &! Array element no.
, lastrow               &! No. of elements in last row
, nlevsrows, nlevscount &! No. of rows and Do Loop counter
, ntroprows, ntropcount &! elements
, nwetrows, nwetcount    ! elements

INTEGER ::              &
  i, j, k, l, m          ! Loop counter

REAL ::                                  &
  alfad0                                 &! Dummy for I/Os
, alfad1(row_length,rows)                &
, alfad2(row_length,rows)                &! Tuning factor for Jan and July
, daynew                                 &! Function to calculate SIN arg
                                          ! (eqn. 2.33 in SCM doc.)
, daysol1,daysol2                        &! No. of days after winter
, dbar0(model_levels)                    &! dummy var for I/Os
, dbar1(row_length,rows,model_levels)    &! Mean dew pt. depressions
, dbar2(row_length,rows,model_levels)    &!  for Jan. and July (K)
, dewpt(row_length,rows,model_levels,2)  &! Dew point (K)
, dgrad0(model_levels)                   &! Dummy var for I/Os
, dgrad1(row_length,rows,model_levels)   &! Gradient dew pt. depressions
, dgrad2(row_length,rows,model_levels)    !  for Jan. and July (K/km)

REAL ::                                 &
  p0(model_levels+1)                    &! Dummy for I/Os
, p1(row_length,rows,model_levels+1)    &! pressure for
, p2(row_length,rows,model_levels+1)    &!  Jan. and July
, p_in(row_length,rows,model_levels+1)  &
                                  ! Initial pressure
, qi(row_length,rows,model_levels)      &! Initial specific humidity (kg/kg)
, rp_theta_levels(row_length,rows,model_levels) &
                                  ! Reciprocal pressure (1/HPa or 1/mb)
, tbar0(model_levels)                   &! dummy var fo I/Os
, tbar1(row_length,rows,model_levels)   &! Mean temp. for Jan. and July
, tbar2(row_length,rows,model_levels)   &!  (K)
, tgrad0(model_levels)                  &! Dummy var for I/Os
, tgrad1(row_length,rows,model_levels)  &! Gradient temp.
, tgrad2(row_length,rows,model_levels)  &!  for Jan. and July (K/km)
, ti(row_length,rows,model_levels)      &! Initial temps. (K)
, tsd0(model_levels)                    &! dummy var fo I/Os
, tsd1(row_length,rows,model_levels)    &! SD of temp. for Jan. and July
, tsd2(row_length,rows,model_levels)    &!  (K)
, tstar0                                &! Dummy for I/Os
, tstar1(row_length,rows)               &! Surface temperature for
, tstar2(row_length,rows)               &!  Jan and July (K)
, vnbar0(model_levels)                  &! Dummy var for I/Os
, vnbar1(row_length,rows,model_levels)  &! Mean horizontal velocity VN
, vnbar2(row_length,rows,model_levels)  &!  for Jan. and July (m/s)
, vnsd0(model_levels)                   &! Dummy for I/Os
, vnsd1(row_length,rows,model_levels)   &! SD horizontal velocity VN
, vnsd2(row_length,rows,model_levels)   &!  for Jan. and July (m/s)
, vpbar0(model_levels)                  &! Dummy var for I/O.
, vpbar1(row_length,rows,model_levels)  &! Mean horizontal velocity VP
, vpbar2(row_length,rows,model_levels)  &!  for Jan. and July (m/s)
, wbar0(ntrop)                          &! Dummy for I/Os
, wbar1(row_length,rows,ntrop)          &! Mean vertical velocity
, wbar2(row_length,rows,ntrop)          &!  for Jan and July (mb/s)
, wsd0(ntrop)                           &! Dummy for I/Os
, wsd1(row_length,rows,ntrop)           &! SD vertical velocity
, wsd2(row_length,rows,ntrop)           &!  for Jan and July (mb/s)
, xt                                     ! Argument of SIN distribution
                                         ! eqn. 2.33

CHARACTER (LEN=10)  :: c1fmt, c2fmt, c3fmt, c4fmt, c5fmt
CHARACTER (LEN=200) :: c6fmt, c7fmt, c8fmt, c9fmt, c10fmt
CHARACTER (LEN=120) :: c11fmt, c12fmt

CHARACTER (LEN=25) :: txt ! A temporary variable used to test
                          ! the format of the file on unit 25
!  Set up Format statements

! The numbers in the files on units 25 and 26 may have been
! written with or without decimal points in them. We need to
! determine which so we can get the format statements right.
! Read the first line of unit 25 into a character variable

READ (25,'(A)') txt
! Does it contain any decimal points?
IF (INDEX(txt,'.') == 0) THEN
  ! No - use the following format statements
  c1fmt = '(6e8.2)'
  c2fmt = '(10f8.2)'
  c3fmt = '(5e8.2)'
  c4fmt = '(4f7.2)'
  c5fmt = '(f7.6)'
ELSE
  ! Yes - assume the numbers are one character longer
  c1fmt = '(6e9.2)'
  c2fmt = '(10f9.2)'
  c3fmt = '(5e9.2)'
  c4fmt = '(4f8.2)'
  c5fmt = '(f8.6)'
END IF

! Rewind unit 25 in readiness for proper reading below
REWIND(25)

c6fmt = '("Climate forcing data for july"'                        &
  //',/," ________________________________"'                      &
  //',/," lat=",f7.2,"  long=",f7.2)'
c7fmt = '("  tstar (K)  ",f7.2/," tuning factor " ,f7.2,/,'       &
  //'" dayno. relative to winter solstice ",f7.2)'
c8fmt = '("0          level",i2,9(4x,"level",i2))'
c9fmt = '('                                                       &
 //'" tmn K    ",10(1e10.3,1x)/," tsd K    ",'                    &
 //'10(1e10.3,1X)/," tgrd K/km",10(1e10.3,1x)/,'                  &
 //'" p mb     ",10(1e10.3,1X)/,'                                 &
 //'" vnmn m/s ",10(1e10.3,1x)/," vpmn m/s ",'                    &
 //'10(1e10.3,1x)/," vnsd m/s ",10(1e10.3,1x))'
c10fmt = '('                                                      &
 //'" dmn K    ",10(1e10.3,1x)/," dgrd K/km",10(1e10.3,1x)/)'
c11fmt = '('                                                      &
 //'" wbar mb/s",10(1e10.3,1x)/," wsd mb/s ",10(1e10.3,1x)/)'
c12fmt = '("1Climate forcing data for january"'                   &
 //',/," ________________________________"'                       &
 //',/," lat=",f7.2,"  long=",f7.2)'

! Read climate stats for January and July
!
! Each column is read in a set of dummy variables,
! then, the values are copied accross to the real
! arrays. This is to ensure conistency between
! one column and multicolumns runs, even though it creates
! redundancy

DO m=1, rows
  DO l=1, row_length
    READ (25,c4fmt) lat0, longitude0, tstar0, alfad0
    READ (25,c2fmt) tbar0
    READ (25,c2fmt) tsd0
    READ (25,c2fmt) p0
    READ (25,c2fmt) dbar0
    READ (25,c3fmt) tgrad0
    READ (25,c3fmt) dgrad0
    READ (25,c2fmt) vnbar0
    READ (25,c2fmt) vpbar0
    READ (25,c2fmt) vnsd0
    READ (25,c1fmt) wbar0
    READ (25,c1fmt) wsd0
    READ (25,c5fmt) daysol1
    lat(l,m) = lat0
    longitude(l,m) = longitude0
    tstar1(l,m) = tstar0
    alfad1(l,m) = alfad0

    DO k=1, model_levels
      tbar1(l,m,k) = tbar0(k)
      tsd1(l,m,k) = tsd0(k)
      p1(l,m,k) = p0(k)
      tgrad1(l,m,k) = tgrad0(k)
      vnbar1(l,m,k) = vnbar0(k)
      vpbar1(l,m,k) = vpbar0(k)
      vnsd1(l,m,k) = vnsd0(k)
    END DO

    p1(l,m,model_levels+1) = p0(model_levels+1)

    DO k=1, model_levels
      dbar1(l,m,k) = dbar0(k)
      dgrad1(l,m,k) = dgrad0(k)
    END DO

    DO k=1, ntrop
      wbar1(l,m,k) = wbar0(k)
      wsd1(l,m,k) = wsd0(k)
    END DO

    READ (26,c4fmt) lat0, longitude0, tstar0, alfad0
    READ (26,c2fmt) tbar0
    READ (26,c2fmt) tsd0
    READ (26,c2fmt) p0
    READ (26,c2fmt) dbar0
    READ (26,c3fmt) tgrad0
    READ (26,c3fmt) dgrad0
    READ (26,c2fmt) vnbar0
    READ (26,c2fmt) vpbar0
    READ (26,c2fmt) vnsd0
    READ (26,c1fmt) wbar0
    READ (26,c1fmt) wsd0
    READ (26,c5fmt) daysol2
    lat(l,m) = lat0
    longitude(l,m) = longitude0
    tstar2(l,m) = tstar0
    alfad2(l,m) = alfad0

    DO k=1, model_levels
      tbar2(l,m,k) = tbar0(k)
      tsd2(l,m,k) = tsd0(k)
      p2(l,m,k) = p0(k)
      tgrad2(l,m,k) = tgrad0(k)
      vnbar2(l,m,k) = vnbar0(k)
      vpbar2(l,m,k) = vpbar0(k)
      vnsd2(l,m,k) = vnsd0(k)
    END DO

    p2(l,m,model_levels+1) = p0(model_levels+1)

    DO k=1, model_levels
      dbar2(l,m,k) = dbar0(k)
      dgrad2(l,m,k) = dgrad0(k)
    END DO

    DO k=1, ntrop
      wbar2(l,m,k) = wbar0(k)
      wsd2(l,m,k) = wsd0(k)
    END DO

  END DO                     ! l
END DO                      ! m

! Calculate amplitude and mean of annual sinusoidal distribution
! Eqs 10 and 11

! DEPENDS ON: abnew
CALL abnew (tstar1, tstar2, tstara, tstarb, row_length, rows, 1)

! DEPENDS ON: abnew
CALL abnew (p1, p2, pa, pb, row_length, rows, model_levels+1)

! DEPENDS ON: abnew
CALL abnew (alfad1, alfad2, alfada, alfadb, row_length, rows, 1)

! DEPENDS ON: abnew
CALL abnew (tbar1, tbar2, tbara, tbarb, row_length, rows, model_levels)

! DEPENDS ON: abnew
CALL abnew (tsd1, tsd2, tsda, tsdb, row_length, rows, model_levels)

! DEPENDS ON: abnew
CALL abnew (tgrad1, tgrad2, tgrada, tgradb, row_length, rows, model_levels)

! DEPENDS ON: abnew
CALL abnew (dbar1, dbar2, dbara, dbarb, row_length, rows, model_levels)

! DEPENDS ON: abnew
CALL abnew (dgrad1, dgrad2, dgrada, dgradb, row_length, rows, model_levels)

! DEPENDS ON: abnew
CALL abnew (vnbar1, vnbar2, vnbara, vnbarb, row_length, rows, model_levels)

! DEPENDS ON: abnew
CALL abnew (vnsd1,  vnsd2,  vnsda,  vnsdb, row_length, rows, model_levels)

! DEPENDS ON: abnew
CALL abnew (vpbar1, vpbar2, vpbara, vpbarb, row_length, rows, model_levels)

! DEPENDS ON: abnew
CALL abnew (wbar1,  wbar2,  wbara,  wbarb, row_length, rows, ntrop)

! DEPENDS ON: abnew
CALL abnew (wsd1,   wsd2,   wsda,   wsdb, row_length, rows, ntrop)

! Calculate constants for annual cycle used in eqn. 12
atime = 2.0 * pi / andayy
btime = pi * (0.5-2.0*daysol1)

! Calculate argument of sinusoidal distribution (eqn. 12)
! DEPENDS ON: daynew
xt = daynew (atime, btime, dayno)

! Calculate sinusoidal distribution (eqn. 12)

! DEPENDS ON: xnew
CALL xnew (p_in, pa, pb, row_length, rows, model_levels+1, xt)

! DEPENDS ON: xnew
CALL xnew(ti, tbara, tbarb, row_length, rows, model_levels, xt)

! DEPENDS ON: xnew
CALL xnew(q, dbara, dbarb, row_length, rows, model_levels, xt)

! Calculate default initial profile for Q

DO k=1, model_levels
  DO j=1, rows
    DO i=1, row_length
      dewpt(i,j,k,1) = ti(i,j,k) - q(i,j,k)
    END DO
  END DO
END DO

CALL qsat                                                                   &
  ( qi, dewpt(:,:,:,1), p_theta_levels,row_length,rows,model_levels)

DO k=1, model_levels
  DO j=1, rows
    DO i=1, row_length
      t(i,j,k) = ti(i,j,k)
    END DO
  END DO
END DO

DO k=1, model_levels
  DO j=1, rows
    DO i=1, row_length
      q(i,j,k) = qi(i,j,k)
    END DO
  END DO
END DO                     ! i

!---------------------------------------------------------------------
!     Print out climate datasets for January and July as read in
!     This section of code is very long but is necessary for
!     flexiblity ie to cope with any number of levels
!---------------------------------------------------------------------

daysol1 = daysol1 * andayy
daysol2 = daysol2 * andayy

DO l=1 , row_length
  DO m=1, rows
    ! Transfer the arrays back to their 1D versions:
    lat0 = lat(l,m)
    longitude0 = longitude(l,m)
    tstar0 = tstar1(l,m)
    alfad0 = alfad1(l,m)

    DO k=1, model_levels
      tbar0(k)  = tbar1(l,m,k)
      tsd0(k)   = tsd1(l,m,k)
      p0(k)     = p1(l,m,k)
      tgrad0(k) = tgrad1(l,m,k)
      vnbar0(k) = vnbar1(l,m,k)
      vpbar0(k) = vpbar1(l,m,k)
      vnsd0(k)  = vnsd1(l,m,k)
    END DO

    DO k=1, model_levels
      dbar0(k)= dbar1(l,m,k)
      dgrad0(k)= dgrad1(l,m,k)
    END DO

    DO k=1, ntrop
      wbar0(k) = wbar1(l,m,k)
      wsd0(k)  = wsd1(l,m,k)
    END DO

    IF (row_length*rows  >   1) THEN
      WRITE (11,*) 'Column no ', l,m
    END IF

    WRITE (11,c6fmt) lat0, longitude0
    WRITE (11,c7fmt) tstar0, alfad0, daysol1

    ! Set format statements
    cfmt = '(''          '',  (1e10.3,1x))'
    ctfmt = '(''0        '',  (3x,''level'',i2,1x))'

    ! Calculate no. of rows and no. of elements in last row for
    ! variables with model_levels

    IF ( MOD(model_levels,10)  ==  0) THEN
      nlevsrows = INT(model_levels/10)
      lastrow = 10
    ELSE
      nlevsrows = INT(model_levels/10) + 1
      lastrow = MOD(model_levels,10)
    END IF

    DO nlevscount = 1, nlevsrows
      element = 10 * (nlevscount-1)

      IF (nlevscount  <   nlevsrows) THEN

        ! Write out all complete rows ie of 10 variables per row

        WRITE (11,c8fmt) (element+i, i = 1, 10)
        WRITE (11,c9fmt) (tbar0(element+i), i = 1, 10),             &
          (tsd0(element+i), i = 1, 10),                             &
          (tgrad0(element+i), i = 1, 10),                           &
          (p0(element+i), i = 1, 10),                               &
          (vnbar0(element+i), i = 1, 10),                           &
          (vpbar0(element+i), i = 1, 10),                           &
          (vnsd0(element+i), i = 1, 10)
      ELSE

        ! Write out last row. Use an internal format statement by
        ! creating a character string. This will enable a variable
        ! format to be created eg NF10.6 where N is the no. of
        ! elements in the last row which can be written into the
        ! format statement via an internal write statement.

        WRITE (ctfmt(14:15),'(i2)')lastrow
        WRITE(11,ctfmt)(element+i,i= 1, lastrow)
        WRITE(cfmt(15:16),'(i2)')lastrow
        WRITE(cfmt(4:12),'(''tmn k    '')')
        WRITE(11,cfmt)(tbar0(i+element),i= 1, lastrow)
        WRITE(cfmt(4:12),'(''tsd k    '')')
        WRITE(11,cfmt)(tsd0(i+element),i= 1, lastrow)
        WRITE(cfmt(4:12),'(''tgrd k/km'')')
        WRITE(11,cfmt)(tgrad0(i+element),i= 1, lastrow)
        WRITE(cfmt(4:12),'(''P mb '')')
        WRITE(11,cfmt)(p0(i+element),i= 1, lastrow)
        WRITE(cfmt(4:12),'(''vnmn m/s '')')
        WRITE(11,cfmt)(vnbar0(i+element),i= 1, lastrow)
        WRITE(cfmt(4:12),'(''vpmn m/s '')')
        WRITE(11,cfmt)(vpbar0(i+element),i= 1, lastrow)
        WRITE(cfmt(4:12),'(''vnsd m/s '')')
        WRITE(11,cfmt)(vnsd0(i+element),i= 1, lastrow)
      END IF
    END DO

    ! Calculate no. of rows and no. of elements in last row for
    ! variables with NWET

    IF ( MOD(model_levels,10)  ==  0) THEN
      nwetrows = INT(model_levels/10)
      lastrow  = 10
    ELSE
      nwetrows = INT(model_levels/10) + 1
      lastrow  = MOD(model_levels,10)
    END IF

    DO nwetcount = 1, nwetrows
      element = 10*(nwetcount-1)
      IF (nwetcount  <   nwetrows) THEN

        ! Write out all complete rows ie of 10 variables per row

        WRITE (11,c8fmt) (element+i,i = 1, 10)
        WRITE (11,c10fmt) (dbar0(element+i),i = 1, 10),             &
          (dgrad0(element+i),i = 1, 10)
      ELSE

        ! Write out last row. Use an internal format statement by
        ! creating a character string. This will enable a variable
        ! format to be created eg NF10.6 where N is the no. of
        ! elements in the last row which can be written into the
        ! format statement via an internal write statement.

        WRITE (ctfmt(14:15),'(i2)') lastrow
        WRITE (11,ctfmt) (element+i, i = 1, lastrow)
        WRITE (cfmt(15:16),'(i2)') lastrow
        WRITE (cfmt(4:12),'(''dmn k    '')')
        WRITE (11,cfmt) (dbar0(i+element), i = 1, lastrow)
        WRITE (cfmt(4:12),'(''dgrd k/km'')')
        WRITE (11,cfmt) (dgrad0(i+element), i = 1, lastrow)
      END IF
    END DO

    ! Calculate no. of rows and no. of elements in last row for
    ! variables with NTROP

    IF ( MOD(ntrop,10)  ==  0) THEN
      ntroprows = INT(ntrop/10)
      lastrow = 10
    ELSE
      ntroprows = INT(ntrop/10) + 1
      lastrow = MOD(ntrop,10)
    END IF

    DO ntropcount = 1,ntroprows
      element = 10 * (ntropcount-1)
      IF ( ntropcount  <   ntroprows) THEN

        ! Write out all complete rows ie of 10 variables per row

        WRITE (11,c8fmt) (element+i,i=1,10)
        WRITE (11,c11fmt) (wbar0(element+i), i = 1, 10),            &
          (wsd0(element+i),i=1,10)
      ELSE

        ! Write out last row. Use an internal format statement by
        ! creating a character string. This will enable a variable
        ! format to be created eg NF10.6 where N is the no. of
        ! elements in the last row which can be written into the
        ! format statement via an internal write statement.

        WRITE (ctfmt(14:15),'(i2)') lastrow
        WRITE (11,ctfmt) (element+i, i = 1, lastrow)
        WRITE (cfmt(15:16),'(i2)') lastrow
        WRITE (cfmt(4:12),'(''wmn mb/s '')')
        WRITE (11,cfmt) (wbar0(i+element), i = 1, lastrow)
        WRITE (cfmt(4:12),'(''wsd mb/s '')')
        WRITE (11,cfmt) (wsd0(i+element), i = 1, lastrow)
      END IF
    END DO

    ! Transfer the arrays back to their 1D versions:
    lat0   = lat(l,m)
    longitude0  = longitude(l,m)
    tstar0 = tstar2(l,m)
    alfad0 = alfad2(l,m)

    DO k=1, model_levels
      tbar0(k)  = tbar2(l,m,k)
      tsd0(k)   = tsd2(l,m,k)
      tgrad0(k) = tgrad2(l,m,k)
      p0(k)     = p2(l,m,k)
      vnbar0(k) = vnbar2(l,m,k)
      vpbar0(k) = vpbar2(l,m,k)
      vnsd0(k)  = vnsd2(l,m,k)
    END DO

    DO k=1, model_levels
      dbar0(k)  = dbar2(l,m,k)
      dgrad0(k) = dgrad2(l,m,k)
    END DO

    DO k=1, ntrop
      wbar0(k) = wbar2(l,m,k)
      wsd0(k)  = wsd2(l,m,k)
    END DO


    WRITE (11,c12fmt) lat0, longitude0
    WRITE (11,c7fmt) tstar0, alfad0, daysol2

    ! Calculate no. of rows and no. of elements in last row for
    ! variables with model_levels

    IF ( MOD(model_levels,10) == 0) THEN
      nlevsrows = INT(model_levels/10)
      lastrow = 10
    ELSE
      nlevsrows = INT(model_levels/10) + 1
      lastrow = MOD(model_levels,10)
    END IF

    DO nlevscount=1, nlevsrows
      element = 10*(nlevscount-1)
      IF (nlevscount  <   nlevsrows) THEN

        ! Write out all complete rows ie of 10 variables per row

        WRITE (11,c8fmt) (element+i, i = 1, 10)
        WRITE (11,c9fmt) (tbar0(element+i),i = 1, 10),              &
          (tsd0(element+i),i = 1, 10),                              &
          (tgrad0(element+i),i = 1, 10),                            &
            (p0(element+i),i = 1, 10),                              &
          (vnbar0(element+i),i = 1, 10),                            &
          (vpbar0(element+i),i = 1, 10),                            &
          (vnsd0(element+i),i = 1, 10)
      ELSE

        ! Write out last row. Use an internal format statement by
        ! a character string. This will enable a variable format
        ! to be created eg NF10.6 where N is the no. of elements in
        ! the last row which can be written into the format
        ! statement via an internal write statement.

        WRITE (ctfmt(14:15),'(i2)') lastrow
        WRITE (11,ctfmt) (element+i,i = 1, lastrow)
        WRITE (cfmt(15:16),'(i2)') lastrow
        WRITE (cfmt(4:12),'(''tmn k    '')')
        WRITE (11,cfmt) (tbar0(i+element), i = 1, lastrow)
        WRITE (cfmt(4:12),'(''tsd k    '')')
        WRITE (11,cfmt) (tsd0(i+element),i = 1, lastrow)
        WRITE (cfmt(4:12),'(''tgrd k/km'')')
        WRITE (11,cfmt) (tgrad0(i+element),i = 1, lastrow)
        WRITE (cfmt(4:12),'(''P mb '')')
        WRITE (11,cfmt) (p0(i+element),i = 1, lastrow)
        WRITE (cfmt(4:12),'(''vnmn m/s '')')
        WRITE (11,cfmt) (vnbar0(i+element),i = 1, lastrow)
        WRITE (cfmt(4:12),'(''vpmn m/s '')')
        WRITE (11,cfmt) (vpbar0(i+element),i = 1, lastrow)
        WRITE (cfmt(4:12),'(''vnsd m/s '')')
        WRITE (11,cfmt) (vnsd0(i+element),i = 1, lastrow)
      END IF
    END DO

    ! Calculate no. of rows and no. of elements in last row for
    ! variables with model_levels

    IF ( MOD(model_levels,10)  ==  0) THEN
      nwetrows = INT(model_levels/10)
      lastrow = 10
    ELSE
      nwetrows = INT(model_levels/10) + 1
      lastrow = MOD(model_levels,10)
    END IF

    DO nwetcount=1, nwetrows
      element = 10 * (nwetcount-1)
      IF (nwetcount  <   nwetrows) THEN

        ! Write out all complete rows ie of 10 variables per row

        WRITE (11,c8fmt) (element+i,i = 1, 10)
        WRITE (11,c10fmt) (dbar0(element+i),i = 1, 10),             &
          (dgrad0(element+i),i = 1, 10)
      ELSE

        ! Write out last row. Use an internal format statement by
        ! creating a character string. This will enable a variable
        ! format to be created eg NF10.6 where N is the no. of
        ! elements in the last row which can be written into the
        ! format statement via an internal write statement.

        WRITE (ctfmt(14:15),'(i2)') lastrow
        WRITE (11,ctfmt) (element+i,i = 1, lastrow)
        WRITE (cfmt(15:16),'(i2)') lastrow
        WRITE (cfmt(4:12),'(''dmn k/km '')')
        WRITE (11,cfmt) (dbar0(i+element),i = 1, lastrow)
        WRITE (cfmt(4:12),'(''dgrd k/km'')')
        WRITE (11,cfmt) (dgrad0(i+element),i = 1, lastrow)
      END IF
    END DO

    ! Calculate no. of rows and no. of elements in last row for
    ! variables with ntrop

    IF ( MOD(ntrop,10) == 0) THEN
      ntroprows = INT(ntrop/10)
      lastrow = 10
    ELSE
      ntroprows = INT(ntrop/10) + 1
      lastrow = MOD(ntrop,10)
    END IF

    DO ntropcount=1, ntroprows
      element = 10*(ntropcount-1)
      IF (ntropcount  <   ntroprows) THEN

        ! Write out all complete rows ie of 10 variables per row

        WRITE (11,c8fmt) (element+i, i = 1, 10)
        WRITE (11,c11fmt) (wbar0(element+i), i = 1, 10),            &
          (wsd0(element+i),i = 1, 10)
      ELSE

        ! Write out last row. Use an internal format statement by
        ! creating a character string. This will enable a variable
        ! formatC to be created eg NF10.6 where N is the no. of
        ! elements in the last row which can be written into the
        ! format statement via an internal write statement.

        WRITE (ctfmt(14:15),'(i2)') lastrow
        WRITE (11,ctfmt)(element+i, i = 1, lastrow)
        WRITE (cfmt(15:16),'(i2)') lastrow
        WRITE (cfmt(4:12),'(''wmn mb/s '')')
        WRITE (11,cfmt)(wbar0(i+element), i = 1, lastrow)
        WRITE (cfmt(4:12),'(''wsd mb/s '')')
        WRITE (11,cfmt)(wsd0(i+element), i = 1, lastrow)
      END IF
    END DO
  END DO                   ! m
END DO                     ! l


RETURN

END SUBROUTINE initstat

