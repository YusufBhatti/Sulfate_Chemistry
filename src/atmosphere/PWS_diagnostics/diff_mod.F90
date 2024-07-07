! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine to calculate diffs wrt p, x and y

MODULE diff_mod

! Description: collection of routines to peform diffs on fields.
!
! Method: splines and finite diffs on input FField
!
! dFdP supports different input grids, ie C or B grid, with suitable inputs.
! dFdx and dFdy only support uv b grid fields.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS diagnostics
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3


 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIFF_MOD'

CONTAINS

SUBROUTINE DiffP( NumLevs,     &  ! in 
                  PRef,         &  ! in
                  ini_start,ini_end,                              &
                  inj_start,inj_end,                              &
                  FFields,      &  ! in
                  PFields,      &  ! in
                  dFdP)            ! out


USE atm_fields_bounds_mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs                      ! No of model levels
REAL, INTENT(IN)    :: PRef                         ! P where deriv reqd

INTEGER, INTENT(IN)  ::   ini_start,ini_end
INTEGER, INTENT(IN)  ::   inj_start,inj_end

REAL, INTENT(IN)    ::  FFields(ini_start:ini_end,                  &
                                inj_start:inj_end,                  &
                                1:NumLevs)  ! Input variable, F

REAL, INTENT(IN)    ::  PFields(ini_start:ini_end,                  &
                                inj_start:inj_end,                  &
                                1:NumLevs) ! Corresp. P vals

REAL, INTENT(OUT)   ::  dfdp(ini_start:ini_end,                              &
                             inj_start:inj_end)    ! Deriv F wrt p

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DIFFP"
INTEGER, PARAMETER :: NPts = 6                      ! No. spline coeffs

! Local Variables:
INTEGER :: i, j, jlwr, jupr, jmid, k, halfn
REAL :: dx, t
REAL :: b(NPts), c(NPts), d(NPts)                       ! Spline coeffs
REAL :: deltay(NPts-1), deltax(NPts-1), dydx(NPts-1)    ! Spline inputs

INTEGER :: Lower (ini_start:ini_end,                              &
                  inj_start:inj_end)    ! Lower model level
INTEGER :: Upper (ini_start:ini_end,                              &
                  inj_start:inj_end)    ! Upper model level

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header --------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL Timer( RoutineName, 3 )

halfn = NPts/2
#if defined(VECTOR)
DO j = inj_start,inj_end
  DO i = ini_start,ini_end
    IF (PRef >= PFields(i,j,halfn)) THEN

       ! all smaller

      Lower(i,j) = halfn
      Upper(i,j) = halfn+1
    ELSE IF (PRef < PFields(i,j,NumLevs+1-halfn)) THEN

       ! all lower

      Lower(i,j) = NumLevs-halfn
      Upper(i,j) =  NumLevs+1-halfn
    ELSE
      Lower (i,j) = -1 ! One level to find
    END IF
  END DO
END DO

DO k=halfn , NumLevs+1-halfn
  DO j = inj_start,inj_end
    DO i = ini_start,ini_end

      IF ( ( PFields(i,j,k) > Pref  ) &
         .AND. ( PFields(i,j,k+1) <= Pref  ) &
         .AND. (Lower(i,j) == -1) ) THEN
        Lower(i,j) = k
        Upper(i,j) = k+1
      END IF

    END DO
  END DO
END DO
#else
DO j = inj_start,inj_end
  DO i =ini_start,ini_end
    ! Divide and conquer to find level containing PRef
    jlwr = halfn
    jupr = NumLevs+1-halfn
    DO WHILE ( (jupr-jlwr) > 1 )
      jmid = (jlwr+jupr)/2
      IF ( PRef >= PFields(i,j,jmid)) THEN
        jupr = jmid
      END IF
      IF ( PRef <  PFields(i,j,jmid)) THEN
        jlwr = jmid
      END IF
    END DO

    Lower (i,j) = jlwr
    Upper (i,j) = jupr

  END DO
END DO
#endif

DO j = inj_start,inj_end
  DO i = ini_start,ini_end
    jlwr = Lower (i,j)
    jupr = Upper (i,j)

    ! Use fields in reverse order so that pressure is increasing
    dx = PRef - PFields(i,j,jupr)
    DO k = 1,NPts-1
      deltax(k) = PFields(i,j,jupr+halfn  -k)-  &
                  PFields(i,j,jupr+halfn-1-k) 
      deltay(k) = FFields(i,j,jupr+halfn  -k)-  &
                  FFields(i,j,jupr+halfn-1-k)
    END DO

    !-----------------------------------------------------------------------
    ! the coefficients b(i), c(i), and d(i), i=1,2,...,NPts are computed
    ! for a cubic interpolating spline
    ! s(u) = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
    !
    !   derivative= b(i) + 2.*c(i)*(u-x(i)) + 3.*d(i)*(u-x(i))**2
    !
    ! where  x(i) < u < x(i+1), using horner's rule
    ! if  u < x(1) then  i = 1  is used.
    ! if  u >= x(n) then  i = n  is used.
    !
    ! A binary search is performed to determine the proper interval.
    ! input..
    !   NPts = the number of data points or knots (NPts >= 2)
    !   x = the abscissas of the knots in strictly increasing order
    !   y = the ordinates of the knots
    ! output..
    !   b, c, d  = arrays of spline coefficients as defined above.
    ! using  p  to denote differentiation,
    !   y(i) = s(x(i))
    !   b(i) = sp(x(i))
    !   c(i) = spp(x(i))/2
    !   d(i) = sppp(x(i))/6  (derivative from the right)
    !-----------------------------------------------------------------------

    ! Tridiagonal system: b = diagonal, d = offdiagonal, c = RH side.
    ! Third derivs at  x(1)  and  x(n) obtained from divided differences
        ! Derivative to the right
    dydx(1:NPts-1) = deltay(1:NPts-1) / deltax(1:NPts-1)
    ! Difference in right and left derivative
    c(2:NPts-1) = dydx  (2:NPts-1) - dydx  (1:NPts-2)
    ! b(k)=x(k+1)-x(k-1)
    b(1)        = -deltax(1)
    b(NPts)     = -deltax(NPts-1)
    b(2:NPts-1) = deltax(1:NPts-2) + deltax(2:NPts-1)

    c(1)    = c(3)     /b(3)      - c(2)     /b(2)
    c(NPts) = c(NPts-2)/b(NPts-2) - c(NPts-1)/b(NPts-1)

    c(1)    = c(1)   *b(1)   **2/( b(3)     -b(1) )
    c(NPts) = c(NPts)*b(NPts)**2/( b(NPts-2)-b(NPts) )

    b(2:NPts-1) = 2.0*b(2:NPts-1)
    !----------------------------------------------  forward elimination
    DO k = 1, NPts-1
      t = deltax(k)/b(k)
      b(k+1) = b(k+1) - t*deltax(k)
      c(k+1) = c(k+1) - t*c(k)
    END DO
    c(NPts) = c(NPts)/b(NPts)
    !------------------------------------------------  back substitution
    DO k = NPts-1, halfn, -1
      c(k) = (c(k) - deltax(k)*c(k+1))/b(k) !2nd derivative / 6
    END DO

    !--------------------------------------------------------------------
    b(halfn) = dydx(halfn) - deltax(halfn)*(c(halfn+1) + 2.0*c(halfn))
    d(halfn) = (c(halfn+1) - c(halfn)) /deltax(halfn) ! = 3rd deriv/6
    c(halfn) = 3.0*c(halfn)                            ! = 2nd deriv/2

    !-------------------------------------------------------------------
    !--- evaluate the derivative of the cubic spline function
    !-----------------------------------------------------------------

    dFdP(i,j) = b(halfn) + dx*(2.0*c(halfn) + 3.0*dx*d(halfn))

  END DO
END DO

CALL Timer( RoutineName, 4 )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE DiffP

! ---------------------------------------------------------------


SUBROUTINE DiffX( FField,               &  ! in
                  dFdX,  fld_type,       &  ! inout
                  ErrorStatus )             ! inout

USE planet_constants_mod, ONLY: planet_radius
USE atm_fields_bounds_mod
USE ereport_mod, ONLY: ereport
USE missing_data_mod, ONLY: rmdi
USE trignometric_mod, ONLY: cos_v_latitude
USE model_domain_mod, ONLY: model_type, mt_lam, l_regular
USE um_parvars
USE um_parparams, ONLY: peast, pwest, pnorth, psouth
USE cderived_mod, ONLY: delta_lambda
USE mpp_conf_mod, ONLY: swap_field_is_vector

USE nlsizes_namelist_mod, ONLY:

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


IMPLICIT NONE

! Subroutine Arguments:

REAL,INTENT(IN)  ::  FField(udims%i_start:udims%i_end,                       &
                            vdims%j_start:vdims%j_end)   ! Input variable, F

REAL,INTENT(OUT) ::  dfdx(udims%i_start:udims%i_end,                         &
                            vdims%j_start:vdims%j_end)   ! Lon derivative of F

INTEGER, INTENT(IN)    :: fld_type
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DIFFX"

! Local Variables:
REAL :: Fhalo(udims_s%i_start:udims_s%i_end,                                &
              vdims_s%j_start:vdims_s%j_end)  ! F with halo

INTEGER :: i,j, i_start, i_end
REAL    :: HdelX                   ! difference weighting factor.
LOGICAL :: rotated

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header --------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Must have a regular, non-rotated lat/lon grid
IF (.NOT. L_regular) THEN
  ErrorStatus=100
  CALL EReport( RoutineName, ErrorStatus,                        &
                "Cannot calculate longitudinal derivative - " // &
                "var res and/or rotated grids are not supported" )
END IF

DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    Fhalo(i,j)=FField(i,j)
  END DO 
END DO

! Call Swapbounds to fill in haloes ready for taking derivatives
CALL swap_bounds(fhalo,udims%i_len,vdims%j_len,1,                &
                 offx,offy,fld_type,swap_field_is_vector)


i_start = udims%i_start
i_end = udims%i_end
dfdx(:,:)=rmdi

IF (model_type == mt_lam) THEN
  IF (at_extremity(PWest)) i_start = udims%i_start+1
  IF (at_extremity(PEast)) i_end = udims%i_end-1
END IF

DO j = vdims%j_start,vdims%j_end  
  DO i = i_start,i_end
    HdelX = 0.5/( Planet_Radius * delta_lambda * &
                    cos_v_latitude(i,j) )
    DFDx(i,j) = HdelX * ( Fhalo(i+1,j) - Fhalo(i-1,j) )
  END DO
END DO


! need to set something sensible at poles (rmdi)
IF (at_extremity(PSouth) ) THEN
  DFDx(:,vdims%j_start) = rmdi
END IF
IF (at_extremity(PNorth) ) THEN
  DFDx(:,vdims%j_end) = rmdi
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE DiffX

! ---------------------------------------------------------------


SUBROUTINE DiffY( FField,                 &  
                   dFdY, fld_type,         &  ! inout
                   ErrorStatus )              ! inout

! Description:  Dfdy
!

USE planet_constants_mod, ONLY: planet_radius
USE atm_fields_bounds_mod

USE model_domain_mod, ONLY: l_regular

USE ereport_mod, ONLY: ereport
USE missing_data_mod, ONLY: rmdi
USE um_parvars
USE um_parparams, ONLY: pnorth, psouth
USE cderived_mod, ONLY: delta_phi
USE mpp_conf_mod, ONLY: swap_field_is_vector

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


IMPLICIT NONE

! Subroutine Arguments:

REAL,INTENT(IN)  ::  FField(udims%i_start:udims%i_end,                         &
                            vdims%j_start:vdims%j_end)     ! Input variable, F

REAL,INTENT(OUT) ::  dfdy(udims%i_start:udims%i_end,                           &
                          vdims%j_start:vdims%j_end)       ! Lat derivative of F

INTEGER, INTENT(IN) :: fld_type
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DIFFY"

! Local Variables:
INTEGER :: i,j, j_start, j_end
REAL :: HdelY  ! difference weighting factor.
LOGICAL :: rotated
REAL :: Fhalo(udims_s%i_start:udims_s%i_end,                   &
              vdims_s%j_start:vdims_s%j_end)  ! F with halo

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle   

! End of header --------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!rotated = (a_realhd(5)/= 90.0 .OR. a_realhd(6) /= 0.0)

! Must have a regular, non-rotated lat/lon grid
IF (.NOT. L_regular) THEN
  ErrorStatus=100
  CALL EReport( RoutineName, ErrorStatus,                        &
                "Cannot calculate latitudinal derivative - " // &
                "var res and/or rotated grids are not supported" )
END IF

DO j = vdims%j_start,vdims%j_end
  DO i = udims%i_start,udims%i_end
    Fhalo(i,j)=FField(i,j)
  END DO 
END DO

! Call Swapbounds to fill in haloes ready for taking derivatives
CALL swap_bounds(fhalo,udims%i_len,vdims%j_len,1,                &
                 offx,offy,fld_type,swap_field_is_vector)

j_start = vdims%j_start
j_end = vdims%j_end
dfdy(:,:)=rmdi

IF (at_extremity(PSouth)) j_start = vdims%j_start+1
IF (at_extremity(PNorth)) j_end = vdims%j_end-1

DO j = j_start,j_end
  DO i = udims%i_start,udims%i_end

    HdelY = 0.5/(Planet_Radius * delta_phi )
    DfDy(i,j) = HdelY * (Fhalo(i,j+1)- Fhalo(i,j-1))
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE DiffY

END MODULE diff_mod
