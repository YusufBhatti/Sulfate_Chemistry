! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE convact_mod

USE, INTRINSIC :: iso_c_binding, ONLY: &
  C_INT64_T, C_DOUBLE, C_F_POINTER, C_LOC, C_BOOL

USE fill_pressure_mod, ONLY: fill_pressure
USE fill_n_dspec_mod, ONLY: fill_n_dspec
USE icaoheight_mod, ONLY: icaoheight

IMPLICIT NONE 

CONTAINS

SUBROUTINE convact( cols,           &
                    rows,           &
                    levels,         &
                    CPNRT_Fld,      &
                    BlkCld_Flds,    &
                    ConCld_Flds,    &
                    Ptheta_Flds,    &
                    rmdi,           &
                    return_heights, &
                    P_CBB_Out,      &
                    P_CBT_Out,      &
                    CBHorE_Out) BIND(C, name="convact")
!
! Description:
!   Calculates CB Action fields
!
! Method:
!  Subroutine ConvAct takes fields of convective cloud combined with the
!  field of convective precipitation rate to produce a mask of Cumulonimbus
!  activity. This Cb mask is then compared to convective cloud top and base
!  pressures to give Cb cloud top and base pressures, and finally converted. 
!  to heights (in Kft) using the ICAO definition.
!
!
!  Contents:
!    P_CBB_Out  = Will contain Cb base pressure or ICAO height
!    P_CBT_Out  = Will contain Cb top pressure or ICAO height
!    CBHorE_Out = Will contain Cb horizontal extent (as index)
!
!  Input :
!    cols           : Number of columns
!    rows           : Number of rows
!    levels         : Number of levels
!    CPNRT_Fld      : Convective Precipitation Rate  (5/205)
!    BlkCld_Flds    : Bulk Cloud Fraction     (0/266)
!    ConCld_Flds    : Convective Cloud Amount (5/212)
!    Ptheta_Flds    : Theta level Pressure    (0/408)
!    rmdi           : Missing data indicator
!    return_heights : If True return ICAO height, False returns pressures
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: WAFC CB Library
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6
!

IMPLICIT NONE

INTEGER(KIND=C_INT64_T), INTENT(IN)  :: cols
INTEGER(KIND=C_INT64_T), INTENT(IN)  :: rows
INTEGER(KIND=C_INT64_T), INTENT(IN)  :: levels
REAL(KIND=C_DOUBLE),     INTENT(IN)  :: CPNRT_Fld(cols, rows)
REAL(KIND=C_DOUBLE),     INTENT(IN)  :: BlkCld_Flds(cols, rows, levels)
REAL(KIND=C_DOUBLE),     INTENT(IN)  :: ConCld_Flds(cols, rows, levels)
REAL(KIND=C_DOUBLE),     INTENT(IN)  :: Ptheta_Flds(cols, rows, levels)
REAL(KIND=C_DOUBLE),     INTENT(IN)  :: rmdi
LOGICAL(KIND=C_BOOL),    INTENT(IN)  :: return_heights

REAL(KIND=C_DOUBLE),     INTENT(OUT), TARGET :: P_CBB_Out(cols, rows)
REAL(KIND=C_DOUBLE),     INTENT(OUT), TARGET :: P_CBT_Out(cols, rows)
REAL(KIND=C_DOUBLE),     INTENT(OUT), TARGET :: CBHorE_Out(cols, rows)

! Intermediery pointers for each of the above
REAL(KIND=C_DOUBLE),  POINTER :: P_CBB_Fld(:, :)
REAL(KIND=C_DOUBLE),  POINTER :: P_CBT_Fld(:, :)
REAL(KIND=C_DOUBLE),  POINTER :: CBHorE_Fld(:, :)

! Work array
REAL :: I_CBT_Fld(cols, rows)
REAL :: I_CBB_Fld(cols, rows)

! Local Variables:
INTEGER :: i, j, k         ! Loop counters
INTEGER :: ii, jj          ! Loop counters in optimisation loops
INTEGER :: z               ! Loop counter for field index
REAL    :: pcnt            ! Holds value when determining Cb horiz extent

! Variables minnm, diffn below used in to check cloud tops
! are above cloud bases
!
REAL :: diffn                   ! difference of CB base pressure
                                !  and CB top pressure
REAL :: minnm                   ! indicates CB base is above CB top if negative

! End of header --------------------------------------------------------

! Associate intermediary pointers
CALL C_F_POINTER (C_LOC(P_CBB_Out), P_CBB_Fld, [cols, rows])
CALL C_F_POINTER (C_LOC(P_CBT_Out), P_CBT_Fld, [cols, rows])
CALL C_F_POINTER (C_LOC(CBHorE_Out), CBHorE_Fld, [cols, rows])

! Initialise output arrays to MDI
P_CBB_Fld(:,:) = rmdi
P_CBT_Fld(:,:) = rmdi
CBHorE_Fld(:,:) = rmdi

I_CBB_Fld(:,:) = rmdi
I_CBT_Fld(:,:) = rmdi

!---- Generate initial Cb and layer cloud masks as follows:
!
! If convective PPTN rate greater than .00025 then Cb. Put into I_CBB_Fld.
DO j = 1, rows
  DO i = 1, cols
    IF ( CPNRT_Fld(i,j)   > 0.00025    ) THEN  
      I_CBB_Fld(i,j) = 1.0
    ELSE
      I_CBB_Fld(i,j) = rmdi
    END IF
  END DO
END DO


!---  Adjust mask fields using subroutine FILL_N_DSPEC to fill in small gaps
!     and smooth edges. Subroutine iterates 30 times.
!
! Note: For Cb mask, need to copy mask field into I_CBT_Fld first
!       as this is used in the iteration process in FILL_N_DSPEC
!       but don't need to do this for the layer cloud mask.

! Cb mask
I_CBT_Fld(:,:) = I_CBB_Fld(:,:)

CALL fill_n_dspec(30, I_CBB_Fld, I_CBT_Fld, rmdi)

I_CBT_Fld(:,:) = I_CBB_Fld(:,:)


!---- Identify pressure at base and top of Cbs ----
! Use:
! I_CBB_Fld  containing Cb mask
! ConCld_Flds(:, :, k)  containing convective cloud fraction at level k
! Ptheta_Flds(:, :, k) containing pressure at level k
!
! Output:
! P_CBB_Fld will contain Cb base pressure
! P_CBT_Fld will contain Cb top pressure
!
! Do this by looping over model levels bottom to top.
!
! For cloud bases:
!  If Cb mask indicates Cb in this cell, and convective cloud amount at this
!  level > 0, and no base has yet been identified, set the Cb base to this
!  level. 
!
! For cloud tops:
!  If Cb mask indicates Cb in this cell, and convective cloud amount at this
!  level > 0, and cloud base has been set and is lower than or at this level
!  then set the cloud top to be this level. This is overwritten each time the
!  model level increases and these conditions are met, until the conditions
!  aren't met anymore.
!

DO k= 1, levels
  DO j = 1, rows
    DO i = 1, cols
      !  find cloud Bases
      IF ( ( ConCld_Flds(i,j,k)  /= 0.0 )  .AND.  &
           ( I_CBB_Fld(i,j)     /= rmdi )  .AND.  &
           ( P_CBB_Fld(i,j)     == rmdi )  ) THEN  
        P_CBB_Fld(i,j) =  Ptheta_Flds(i,j,k)
      END IF

      ! Find cloud tops
      IF ( ( ConCld_Flds(i,j,k)  /= 0.0 ) .AND.  &
           ( I_CBB_Fld(i,j)     /= rmdi ) .AND.  &
           ( P_CBB_Fld(i,j) >= Ptheta_Flds(i,j,k))) THEN
        P_CBT_Fld(i,j) = Ptheta_Flds(i,j,k)
      END IF
    END DO
  END DO
END DO



!---- Test that Cb tops are above bases ----
!
! Do this by calculating the difference in pressure between the cloud base
!  and top. If this is negative, the cloud base is above cloud top. In this
!  case set the flag minnm to be equal to the difference. This is then used to 
!  reset all values in the Cb base pressure to the min value. (I.e. if any one 
!  point is wrong, all values in the field are set to the min value.)

minnm=1.0

DO j = 1, rows
  DO i = 1, cols
    diffn = P_CBB_Fld(i,j) - P_CBT_Fld(i,j)
    IF (diffn < 0.0) THEN
      minnm = diffn
    END IF
  END DO
END DO

DO j = 1, rows
  DO i = 1, cols
    IF (minnm < 0.0) THEN
      P_CBB_Fld(i,j) = minnm
    END IF
  END DO
END DO


!---- Apply constraints on Cb top height ----
! Constrain Cb tops such that Cb top must be higher than 400hPa.
! Removes excessively low tops and brings in line with WAFC Washington methods.
! If CB top is lower than 400hPa (i.e. pressure at CB top is > 40000.0Pa )
! then cloud is not assumed  to be Cb and is replaced with RMDI.
! CB mask is also set to RMDI.
DO j = 1, rows
  DO i = 1, cols
    IF ( ( P_CBT_Fld(i,j) > 40000.0 ) .AND.  &
         ( I_CBB_Fld(i,j)/= rmdi ) ) THEN
      P_CBT_Fld(i,j) = rmdi
      P_CBB_Fld(i,j) = rmdi
      I_CBB_Fld(i,j) = rmdi
    END IF
  END DO
END DO


!---- Apply constraints on Cb depth ----
! Ensure Cb depth is realistic, must be greater than or equal to 300hPa.
! If cloud depth not reached, cloud is not assumed  to be Cb and is
! replaced with RMDI. CB mask is also set to RMDI.
DO j = 1, rows
  DO i = 1, cols
    IF (( ( P_CBB_Fld(i,j) ) - ( P_CBT_Fld(i,j) ) < 30000.0 ) &
          .AND. ( P_CBT_Fld(i,j) /= rmdi ) .AND. &
          ( P_CBB_Fld(i,j) /= rmdi )) THEN
      P_CBT_Fld(i,j) = rmdi
      P_CBB_Fld(i,j) = rmdi
      I_CBB_Fld(i,j) = rmdi
    END IF
  END DO
END DO


!---- Calculate Cb horizontal extent ----
! Input:
!   I_CBB_Fld    contains Cb mask
!   ConCld_Flds(:,:,k) contains convective cloud fraction at level k
!
! Output:
!   CBHorE_Fld    containing Cb horizontal extent, given as value of convective
!                 cloud fraction.
!
!  For each gridpoint:
!  1) Set CBHorE_Fld=largest convective cloud fraction in the column
!

DO k = 1, levels
  DO j = 1, rows
    DO i = 1, cols
      IF ( (ConCld_Flds(i,j,k) > CBHorE_Fld(i,j)) .AND. &
           (I_CBB_Fld(i,j) > 0.0 )  ) THEN
        CBHorE_Fld(i,j) = ConCld_Flds(i,j,k)
      END IF
    END DO
  END DO
END DO

DO j = 1, rows
  DO i = 1, cols
    pcnt = rmdi
    IF ( CBHorE_Fld(i,j) > 0.0 ) THEN
      pcnt = CBHorE_Fld(i,j)
    END IF
    CBHorE_Fld(i,j) = pcnt
  END DO
END DO

!---- Fill in field values ----
!
!
! Inputs/output for subroutine Fill_Pressure:
! In_field  - index of input field to be filled in. This will be adjusted
!               within subroutine. On output will hold the adjusted field.
! In_mask   - index of mask field used for comparison (will be either 3 or 7)
! Out_field - index of dummy field which will be overwritten and used within
!               subroutine to compare field values.
!
! On output, the adjusted field will be in both In_field and
! Out_field, and either could be used. The code below uses
! In_field for the adjusted field and Out_field as a
! dummy field which is overwritten.
!

!   Fill in Cloud Base Pressure values within area of Cb mask.
! P_CBB_Fld : Field to be adjusted (will be replaced by adjusted field)
! I_CBB_Fld : Mask field
! I_CBT_Fld : Dummy field, will be overwritten in the subroutine
!             and replaced by adjusted field

CALL fill_pressure(P_CBB_Fld, I_CBB_Fld, I_CBT_Fld, rmdi)


!   Fill in Cloud Top Pressure values within area of Cb mask.
! P_CBT_Fld : Field to be adjusted (will be replaced by adjusted field)
! I_CBB_Fld : Mask field
! I_CBT_Fld : Dummy field, will be overwritten in the subroutine
!             and replaced by adjusted field

CALL fill_pressure(P_CBT_Fld, I_CBB_Fld, I_CBT_Fld, rmdi)

!   Fill in Cb horizontal extent values within area of Cb mask.
! CBHorE_Fld : Field to be adjusted (will be replaced by adjusted field)
! I_CBB_Fld  : Mask field
! I_CBT_Fld  : Dummy field, will be overwritten in the subroutine
!              and replaced by adjusted field

CALL fill_pressure(CBHorE_Fld, I_CBB_Fld, I_CBT_Fld, rmdi)


! Convert the pressures to ICAO heights if requested
IF (return_heights) THEN
  CALL icaoheight(P_CBB_Fld, rmdi)
  CALL icaoheight(P_CBT_Fld, rmdi)
END IF

! Tidy up pointers
NULLIFY(P_CBB_Fld)
NULLIFY(P_CBT_Fld)
NULLIFY(CBHorE_Fld)

END SUBROUTINE convact

END MODULE convact_mod
