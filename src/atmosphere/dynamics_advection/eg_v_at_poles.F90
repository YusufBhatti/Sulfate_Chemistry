! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_v_at_poles_mod

#if !defined(RECON)

USE global_2d_sums_mod, ONLY: &
    global_2d_sums

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_V_AT_POLES_MOD'

CONTAINS
SUBROUTINE eg_v_at_poles(u,v, pole_sgn_flip,                                 &
                   jj_u_pole, jj_v_pole, udim_in,vdim_in,u_pole_v,u_pole_u,  &
                   v_pole_out,xi1_pole_out)


USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE conversions_mod, ONLY: pi
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE um_parvars,  ONLY: gc_proc_row_group
USE eg_parameters_mod, ONLY: pole_consts
USE ereport_mod, ONLY: ereport

IMPLICIT NONE
!
! Description: Calcualte v at the north and south poles
!
!
! Method: Chapter 8, ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_V_AT_POLES'


REAL,    INTENT(IN) :: pole_sgn_flip
INTEGER, INTENT(IN) :: jj_v_pole, jj_u_pole

TYPE (array_dims)   :: udim_in,vdim_in

REAL, INTENT(INOUT) :: v(vdim_in%i_start:vdim_in%i_end,                        &
                         vdim_in%j_start:vdim_in%j_end,                        &
                         vdim_in%k_start:vdim_in%k_end)

REAL, INTENT(IN)    :: u(udim_in%i_start:udim_in%i_end,                        &
                         udim_in%j_start:udim_in%j_end,                        &
                         udim_in%k_start:udim_in%k_end)

! u-component on pole on v-staggered (or p-like) point
REAL, OPTIONAL, INTENT(INOUT)  :: u_pole_v(vdims%i_start:vdims%i_end,1,        &
                                           vdims%k_start:vdims%k_end)

! u-component at pole on u-staggered point
REAL, OPTIONAL, INTENT(INOUT)  :: u_pole_u(udims%i_start:udims%i_end,1,        &
                                           udims%k_start:udims%k_end)

REAL, OPTIONAL, INTENT(OUT) ::   v_pole_out(udim_in%k_start:udim_in%k_end)
REAL, OPTIONAL, INTENT(OUT) :: xi1_pole_out(udim_in%k_start:udim_in%k_end)

! Local variables

INTEGER :: i, k, info
REAL    :: c, d, e, f, dx, pi_fac
REAL    :: a, b, g, xi1_pole, v_pole
REAL    :: glob_sum(udims%k_start:udims%k_end,3)
REAL    :: temp1, temp2, cs, sn

#if defined (INTEL_FORTRAN) && (INTEL_FORTRAN < 13001003)
! There is a bug with Intel Fortran 12.0 which means the direct computation
! fails. The alternative needs a differently shaped array. #5044
REAL    :: tmp_loc(udims%i_start:udims%i_end,udims%k_start:udims%k_end)
#else
REAL    :: tmp_loc(udims%i_start:udims%i_end, udims%k_start:udims%k_end,3)
#endif

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

pi_fac = 1.0/pi

! ************************************************************
! START OF SHARED CODE SECTION (with RCF - below)
!
! First do south pole

c = pole_consts(1)
d = pole_consts(2)
e = pole_consts(3)
f = pole_consts(4)
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( udims, pi_fac, xi1_p, tmp_loc, u, jj_u_pole, xi1_u )     &
!$OMP PRIVATE( i, k, dx )
DO k = udims%k_start, udims%k_end
  DO i = udims%i_start, udims%i_end
    dx = pi_fac*( xi1_p(i+1) - xi1_p(i) )
#if defined (INTEL_FORTRAN) && (INTEL_FORTRAN < 13001003)
    ! There is a bug with Intel Fortran 12.0 which means the direct
    ! computation fails. The alternative needs a differently shaped array.
    ! #5044
    tmp_loc(i,k) = dx*u(i,jj_u_pole,k)*SIN(xi1_u(i))
#else
    tmp_loc(i,k,1) = dx*u(i,jj_u_pole,k)*SIN(xi1_u(i))
    tmp_loc(i,k,2) = dx*u(i,jj_u_pole,k)*COS(xi1_u(i))
    tmp_loc(i,k,3) = dx*u(i,jj_u_pole,k)
#endif
  END DO
END DO
!$OMP END PARALLEL DO

#if defined (INTEL_FORTRAN) && (INTEL_FORTRAN < 13001003)
! There is a bug with Intel Fortran 12.0 which means the direct computation
! fails. The alternative needs a differently shaped array. #5044
CALL global_2d_sums(tmp_loc(:,:),udims%i_len,        &
                    1, 0, 0,  udims%k_len ,          &
                    glob_sum(:,1),  gc_proc_row_group)

DO k = udims%k_start, udims%k_end
  DO i = udims%i_start, udims%i_end
    dx = pi_fac*( xi1_p(i+1) - xi1_p(i) )
    tmp_loc(i,k) = dx*u(i,jj_u_pole,k)*COS(xi1_u(i))
  END DO
END DO



CALL global_2d_sums(tmp_loc(:,:),udims%i_len,        &
                    1, 0, 0,  udims%k_len ,          &
                    glob_sum(:,2),  gc_proc_row_group)

DO k = udims%k_start, udims%k_end
  DO i = udims%i_start, udims%i_end
    dx = pi_fac*( xi1_p(i+1) - xi1_p(i) )
    tmp_loc(i,k) = dx*u(i,jj_u_pole,k)
  END DO
END DO

CALL global_2d_sums(tmp_loc(:,:),udims%i_len,        &
                    1, 0, 0,  udims%k_len ,          &
                    glob_sum(:,3),  gc_proc_row_group)
#else
CALL global_2d_sums(tmp_loc,udims%i_len,             &
                    1, 0, 0, 3*udims%k_len,          &
                    glob_sum,  gc_proc_row_group)
#endif

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( udims, glob_sum, c, d, e, f, pole_sgn_flip, xi2_p,       &
!$OMP         jj_u_pole, vdims, v, jj_v_pole, xi1_p, u_pole_v,         &
!$OMP         u_pole_u, xi1_u, v_pole_out, xi1_pole_out )              &
!$OMP PRIVATE( i, k, a, b, g, temp1, temp2, xi1_pole, cs, sn, dx,      &
!$OMP          v_pole )
DO k = udims%k_start, udims%k_end
  a     = glob_sum(k,1)
  b     = glob_sum(k,2)
  g     = glob_sum(k,3)

  temp1 = (1.0-c-2.0*e**2)*(b-f*g) - (a-e*g)*(d-2.0*e*f)
  temp2 =-(1.0+c-2.0*f**2)*(a-e*g) + (b-f*g)*(d-2.0*e*f)

  temp1 = pole_sgn_flip*temp1
  temp2 = pole_sgn_flip*temp2

  IF ( temp1 == 0.0 .AND. temp2 == 0.0 ) THEN
    xi1_pole = 0.0
  ELSE
    xi1_pole = ATAN2(temp1,temp2)
  END IF

  cs     = COS(xi1_pole)
  sn     = SIN(xi1_pole)
  dx     = SIN(xi2_p(jj_u_pole))
  temp1  = dx*(1.0 - c*(1.0-2.0*sn**2) - 2.0*d*sn*cs                 &
                   - 2.0*(e*cs - f*sn)**2 )
  v_pole = ( (a-e*g)*cs - (b-f*g)*sn )/temp1

  DO i = vdims%i_start, vdims%i_end
    v(i,jj_v_pole,k) = v_pole*COS(xi1_p(i)-xi1_pole)
  END DO

  !
  ! END OF SHARED CODE SECTION
  ! ***************************************************************************

  IF (PRESENT(u_pole_v)) THEN
    DO i = vdims%i_start, vdims%i_end
      u_pole_v(i,1,k) =  -1.0*pole_sgn_flip*v_pole*SIN(xi1_p(i)-xi1_pole)
    END DO

  END IF

  IF (PRESENT(u_pole_u)) THEN

    DO i = udims%i_start, udims%i_end
      u_pole_u(i,1,k) = -1.0*pole_sgn_flip*v_pole*SIN(xi1_u(i)-xi1_pole)
    END DO

  END IF

  IF (PRESENT(  v_pole_out))   v_pole_out(k)=   v_pole
  IF (PRESENT(xi1_pole_out)) xi1_pole_out(k)= xi1_pole


END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_v_at_poles


! =====================================================================

SUBROUTINE eg_uv_at_poles_Bgrid &
                  (u, row_length, rows, n_rows,                       &
                   offx, offy, pole_consts, pole_sgn_flip,            &
                   j_pole_src, j_pole_dest, v)


USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE conversions_mod, ONLY: pi
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE um_parvars, ONLY:  gc_proc_row_group

IMPLICIT NONE
!
! Description: Calculate u or v at the north and south poles on
!              B grid pressure levels, one level at a time
!
!
! Method: Chapter 8, ENDGame formulation version 1.01
!
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_UV_AT_POLES_BGRID'


REAL,    INTENT(IN) :: pole_sgn_flip
INTEGER, INTENT(IN) :: j_pole_src, j_pole_dest

! Array dimensions
INTEGER :: offx, offy
INTEGER :: row_length, rows, n_rows

REAL :: pole_consts(4)

! u and v on B grid
! Note: if v present, v is computed; otherwise u is computed
REAL, INTENT(INOUT) :: u(udims%i_start:udims%i_end, vdims%j_start:vdims%j_end)
REAL, INTENT(INOUT), OPTIONAL :: &
                       v(udims%i_start:udims%i_end, vdims%j_start:vdims%j_end)

! Local variables

INTEGER :: i, info
REAL    :: c, d, e, f, dx, pi_fac
REAL    :: a, b, g, xi1_pole
REAL    :: mag_pole_wind  ! Magnitude of polar wind
REAL    :: glob_sum(3)
REAL    :: temp1, temp2, cs, sn
REAL    :: tmp_loc(udims%i_start:udims%i_end,3)

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

pi_fac = 1.0/pi

c = pole_consts(1)
d = pole_consts(2)
e = pole_consts(3)
f = pole_consts(4)

tmp_loc(:,:) = 0.0

DO i = udims%i_start, udims%i_end
  dx = pi_fac*( xi1_p(i+1) - xi1_p(i) )
  tmp_loc(i,1) = dx *u(i,j_pole_src) *SIN(xi1_u(i))
  tmp_loc(i,2) = dx *u(i,j_pole_src) *COS(xi1_u(i))
  tmp_loc(i,3) = dx *u(i,j_pole_src)
END DO


CALL global_2d_sums(tmp_loc,row_length,        &
                    1, 0, 0, 3,                &
                    glob_sum,  gc_proc_row_group)

a     = glob_sum(1)
b     = glob_sum(2)
g     = glob_sum(3)

temp1 = (1.0-c-2.0*e**2)*(b-f*g) - (a-e*g)*(d-2.0*e*f)
temp2 =-(1.0+c-2.0*f**2)*(a-e*g) + (b-f*g)*(d-2.0*e*f)

temp1 = pole_sgn_flip*temp1
temp2 = pole_sgn_flip*temp2

IF ( temp1 == 0.0 .AND. temp2 == 0.0 ) THEN
  xi1_pole = 0.0
ELSE
  xi1_pole = ATAN2(temp1,temp2)
END IF

cs     = COS(xi1_pole)
sn     = SIN(xi1_pole)
dx     = SIN(xi2_p(j_pole_src))
temp1  = dx*(1.0 - c*(1.0-2.0*sn**2) - 2.0*d*sn*cs                 &
                 - 2.0*(e*cs - f*sn)**2 )
mag_pole_wind = ( (a-e*g)*cs - (b-f*g)*sn )/temp1

DO i = udims%i_start, udims%i_end
  IF (PRESENT(v)) THEN
    v(i,j_pole_dest) = mag_pole_wind*COS(xi1_u(i)-xi1_pole)
  ELSE
    u(i,j_pole_dest) = pole_sgn_flip*mag_pole_wind*SIN(xi1_u(i)-xi1_pole)
  END IF
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_uv_at_poles_Bgrid

#else
IMPLICIT NONE


TYPE array_dims
  INTEGER :: i_start =-HUGE(INT(1))
  INTEGER :: i_end   = HUGE(INT(1))
  INTEGER :: j_start =-HUGE(INT(1))
  INTEGER :: j_end   = HUGE(INT(1))
  INTEGER :: k_start =-HUGE(INT(1))
  INTEGER :: k_end   = HUGE(INT(1))
  INTEGER :: halo_i  = HUGE(INT(1))
  INTEGER :: halo_j  = HUGE(INT(1))
END TYPE array_dims

CONTAINS
SUBROUTINE eg_v_at_poles(u,v,xi1_u, xi1_p, xi2_p, pole_sgn_flip,      &
                             jj_u_pole, jj_v_pole, udim_in,vdim_in)
USE conversions_mod, ONLY: pi


IMPLICIT NONE


! Description:
!
!  an implementation of eg_v_at_poles.F90

! Method:

! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


REAL,    INTENT(IN) :: pole_sgn_flip
INTEGER, INTENT(IN) :: jj_v_pole, jj_u_pole

TYPE (array_dims)   :: udim_in,vdim_in

REAL, INTENT(INOUT) :: v(vdim_in%i_start:vdim_in%i_end,                        &
                         vdim_in%j_start:vdim_in%j_end,                        &
                         vdim_in%k_start:vdim_in%k_end)

REAL, INTENT(IN)    :: u(udim_in%i_start:udim_in%i_end,                        &
                         udim_in%j_start:udim_in%j_end,                        &
                         udim_in%k_start:udim_in%k_end)

REAL, INTENT(IN) :: xi1_u(udim_in%i_start:udim_in%i_end)
REAL, INTENT(IN) :: xi1_p(vdim_in%i_start:vdim_in%i_end)
REAL, INTENT(IN) :: xi2_p(udim_in%j_start:udim_in%j_end)


! Local variables

INTEGER :: i, k, info
REAL    :: c, d, e, f, dx, pi_fac
REAL    :: a, b, g, xi1_pole, v_pole
REAL    :: glob_sum(udim_in%k_start:udim_in%k_end,3)
REAL    :: temp1, temp2, cs, sn
REAL    :: tmp_loc(udim_in%i_start:udim_in%i_end, &
                   udim_in%k_start:udim_in%k_end,3)

REAL :: pole_consts(4)


pi_fac = 1.0/pi

i = udim_in%i_start
dx = ( xi1_p(i+1) - xi1_p(vdim_in%i_end) + 2.0*pi)/(2.0*pi)
pole_consts (1) = dx*COS(2.0*xi1_u(i))
pole_consts (2) = dx*SIN(2.0*xi1_u(i))
pole_consts (3) = dx*SIN(    xi1_u(i))
pole_consts (4) = dx*COS(    xi1_u(i))

DO i = udim_in%i_start+1, udim_in%i_end
  dx = ( xi1_p(i+1) - xi1_p(i) )/(2.0*pi)
  pole_consts (1) = pole_consts (1) + dx*COS(2.0*xi1_u(i))
  pole_consts (2) = pole_consts (2) + dx*SIN(2.0*xi1_u(i))
  pole_consts (3) = pole_consts (3) + dx*SIN(    xi1_u(i))
  pole_consts (4) = pole_consts (4) + dx*COS(    xi1_u(i))
END DO


! First do south pole

c = pole_consts(1)
d = pole_consts(2)
e = pole_consts(3)
f = pole_consts(4)
DO k = udim_in%k_start, udim_in%k_end

  i = udim_in%i_start
  dx = pi_fac*( xi1_p(i+1) -  xi1_p(vdim_in%i_end) + 2.0*pi )
  tmp_loc(i,k,1) = dx*u(i,jj_u_pole,k)*SIN(xi1_u(i))
  tmp_loc(i,k,2) = dx*u(i,jj_u_pole,k)*COS(xi1_u(i))
  tmp_loc(i,k,3) = dx*u(i,jj_u_pole,k)

  DO i = udim_in%i_start+1, udim_in%i_end
    dx = pi_fac*( xi1_p(i+1) - xi1_p(i) )
    tmp_loc(i,k,1) = dx*u(i,jj_u_pole,k)*SIN(xi1_u(i))
    tmp_loc(i,k,2) = dx*u(i,jj_u_pole,k)*COS(xi1_u(i))
    tmp_loc(i,k,3) = dx*u(i,jj_u_pole,k)
  END DO
END DO


glob_sum(udim_in%k_start:udim_in%k_end,1) = &
     SUM(tmp_loc(:,udim_in%k_start:udim_in%k_end,1) )
glob_sum(udim_in%k_start:udim_in%k_end,2) = &
     SUM(tmp_loc(:,udim_in%k_start:udim_in%k_end,2) )
glob_sum(udim_in%k_start:udim_in%k_end,3) = &
     SUM(tmp_loc(:,udim_in%k_start:udim_in%k_end,3) )


DO k = udim_in%k_start, udim_in%k_end
  a     = glob_sum(k,1)
  b     = glob_sum(k,2)
  g     = glob_sum(k,3)

  temp1 = (1.0-c-2.0*e**2)*(b-f*g) - (a-e*g)*(d-2.0*e*f)
  temp2 =-(1.0+c-2.0*f**2)*(a-e*g) + (b-f*g)*(d-2.0*e*f)

  temp1 = pole_sgn_flip*temp1
  temp2 = pole_sgn_flip*temp2

  IF ( temp1 == 0.0 .AND. temp2 == 0.0 ) THEN
    xi1_pole = 0.0
  ELSE
    xi1_pole = ATAN2(temp1,temp2)
  END IF

  cs     = COS(xi1_pole)
  sn     = SIN(xi1_pole)
  dx     = SIN(xi2_p(jj_u_pole))
  temp1  = dx*(1.0 - c*(1.0-2.0*sn**2) - 2.0*d*sn*cs                 &
                   - 2.0*(e*cs - f*sn)**2 )
  v_pole = ( (a-e*g)*cs - (b-f*g)*sn )/temp1

  DO i = vdim_in%i_start, vdim_in%i_end
    v(i,jj_v_pole,k) = v_pole*COS(xi1_p(i)-xi1_pole)
  END DO


END DO

RETURN
END SUBROUTINE eg_v_at_poles
#endif

END MODULE eg_v_at_poles_mod
