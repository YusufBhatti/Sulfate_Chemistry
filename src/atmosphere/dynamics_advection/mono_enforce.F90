! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Mono_Enforce.

SUBROUTINE Mono_Enforce(                                          &
                        Ext_Data, number_of_inputs,               &
                        dim_i_in, dim_j_in, dim_k_in,             &
                        dim_i_out, dim_j_out, dim_k_out,          &
                        halo_i, halo_j,                           &
                        i_out, j_out, k_out,                      &
                        Data_out_high)

! Purpose:
!          Ensures scheme is monotone.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE um_types, ONLY: integer32
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

INTEGER ::                                                        &
  dim_i_in,                                                       &
              ! Dimension of input Data arrays in i direction.
  dim_j_in,                                                       &
              ! Dimension of input Data arrays in j direction.
  dim_k_in,                                                       &
              ! Dimension of input Data arrays in k direction.
  dim_i_out,                                                      &
              ! Dimension of output Data arrays in i direction.
  dim_j_out,                                                      &
              ! Dimension of output Data arrays in j direction.
  dim_k_out,                                                      &
              ! Dimension of output Data arrays in k direction.
  halo_i,                                                         &
              ! Size of halo in i direction.
  halo_j,                                                         &
              ! Size of halo in j direction.
  number_of_inputs

INTEGER (KIND=integer32) ::                                       &
  i_out (dim_i_out* dim_j_out* dim_k_out),                        &
  j_out (dim_i_out* dim_j_out* dim_k_out),                        &
  k_out (dim_i_out* dim_j_out* dim_k_out)

REAL ::                                                           &
  Ext_Data (1-halo_i:dim_i_in+halo_i+1,                           &
            1-halo_j:dim_j_in+halo_j,-1:dim_k_in+2                &
            ,number_of_inputs) !data to be interpolated


! Arguments with Intent IN/OUT.

REAL ::                                                           &
  Data_out_high (dim_i_out* dim_j_out* dim_k_out,                 &
                 number_of_inputs)

! Local Variables.

! scalars

INTEGER ::                                                        &
  ind, n        ! Loop indices

INTEGER :: dim_out

REAL ::                                                           &
  max_mono,                                                       &
  min_mono

INTEGER :: ii      ! Indexing
INTEGER :: jj      ! Indexing
INTEGER :: kk      ! Indexing
REAL :: t1         ! Temporary
REAL :: t2         ! Temporary
REAL :: t3         ! Temporary
REAL :: t4         ! Temporary
REAL :: t5         ! Temporary
REAL :: t6         ! Temporary
REAL :: t7         ! Temporary
REAL :: t8         ! Temporary
REAL :: tmax1      ! Temporary max
REAL :: tmax2      ! Temporary max
REAL :: tmax3      ! Temporary max
REAL :: tmax4      ! Temporary max
REAL :: tmin1      ! Temporary min
REAL :: tmin2      ! Temporary min
REAL :: tmin3      ! Temporary min
REAL :: tmin4      ! Temporary min

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MONO_ENFORCE'


! External Routines: None
! subroutines: None
! Functions: None

! ----------------------------------------------------------------------
!  Section 1.  Set alpha_max to 1 as both input fields are monotone.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

dim_out = dim_i_out*dim_j_out*dim_k_out

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(n,ind,max_mono)                    &
!$OMP& PRIVATE(min_mono, ii,jj,kk, t1,t2,t3,t4,t5,t6,t7,t8)             &
!$OMP& PRIVATE(tmax1,tmax2,tmax3,tmax4, tmin1,tmin2,tmin3,tmin4)        &
!$OMP& SHARED(number_of_inputs,dim_out, Ext_Data,                       &
!$OMP&  data_out_high,i_out,j_out,k_out) 
!DIR$ NOBLOCKING
DO n = 1, number_of_inputs
!$OMP DO SCHEDULE(STATIC) 
#if defined(CRAY_FORTRAN)
!DIR$ VECTOR ALWAYS
!DIR$ SAFE_CONDITIONAL
#endif
  DO ind = 1, dim_out

        ! Setup addressing
        ii = i_out(ind)
        jj = j_out(ind)
        kk = k_out(ind)

        ! Find max and min monotone values for the point concerned.

#if defined(IBM_XL_FORTRAN)
        ! Version optimised for IBM by using the fsel IBM-only intrinsic

        ! To simplify later expressions assign data to temporaries
        t1 = Ext_Data(ii  ,jj  ,kk  ,n)
        t2 = Ext_Data(ii+1,jj  ,kk  ,n)
        t3 = Ext_Data(ii  ,jj+1,kk  ,n)
        t4 = Ext_Data(ii+1,jj+1,kk  ,n)
        t5 = Ext_Data(ii  ,jj  ,kk+1,n)
        t6 = Ext_Data(ii+1,jj  ,kk+1,n)
        t7 = Ext_Data(ii  ,jj+1,kk+1,n)
        t8 = Ext_Data(ii+1,jj+1,kk+1,n)

        ! Max and min are found by a heirarchy of comparisons
        ! Top level comparisons
        tmax1 = fsel(t1-t2,t1,t2)
        tmin1 = fsel(t1-t2,t2,t1)
        tmax2 = fsel(t3-t4,t3,t4)
        tmin2 = fsel(t3-t4,t4,t3)
        tmax3 = fsel(t5-t6,t5,t6)
        tmin3 = fsel(t5-t6,t6,t5)
        tmax4 = fsel(t7-t8,t7,t8)
        tmin4 = fsel(t7-t8,t8,t7)

        ! Second level comparisons
        tmax1 = fsel(tmax1-tmax2,tmax1,tmax2)
        tmin1 = fsel(tmin1-tmin2,tmin2,tmin1)
        tmax2 = fsel(tmax3-tmax4,tmax3,tmax4)
        tmin2 = fsel(tmin3-tmin4,tmin4,tmin3)

        ! Final comparisons
        tmax1 = fsel(tmax1-tmax2,tmax1,tmax2)
        tmin1 = fsel(tmin1-tmin2,tmin2,tmin1)

        max_mono = tmax1
        min_mono = tmin1

        IF ( data_out_high(ind,n)  >   max_mono )               &
             data_out_high(ind,n) = max_mono
        IF ( data_out_high(ind,n)  <   min_mono )               &
             data_out_high(ind,n) = min_mono

#elif  defined(CRAY_FORTRAN)

        max_mono = MAX (                                         &
          Ext_Data(ii,jj,  kk,  n), Ext_Data(ii+1,jj,  kk,  n),  &
          Ext_Data(ii,jj+1,kk,  n), Ext_Data(ii+1,jj+1,kk,  n),  &
          Ext_Data(ii,jj,  kk+1,n), Ext_Data(ii+1,jj,  kk+1,n),  &
          Ext_Data(ii,jj+1,kk+1,n), Ext_Data(ii+1,jj+1,kk+1,n) )

        min_mono = MIN (                                         &
          Ext_Data(ii,jj,  kk,  n), Ext_Data(ii+1,jj,  kk,  n),  &
          Ext_Data(ii,jj+1,kk,  n), Ext_Data(ii+1,jj+1,kk,  n),  &
          Ext_Data(ii,jj,  kk+1,n), Ext_Data(ii+1,jj,  kk+1,n),  &
          Ext_Data(ii,jj+1,kk+1,n), Ext_Data(ii+1,jj+1,kk+1,n) )

        IF (data_out_high(ind,n) > Ext_Data(ii,jj,  kk,  n)) THEN

          IF ( data_out_high(ind,n)  >   max_mono )            &
               data_out_high(ind,n) = max_mono
        ELSE

          IF ( data_out_high(ind,n)  <   min_mono )           &
               data_out_high(ind,n) = min_mono
        END IF

#else

        IF (data_out_high(ind,n) > Ext_Data(ii,jj,  kk,  n)) THEN

          max_mono = MAX (                                       &
          Ext_Data(ii,jj,  kk,  n), Ext_Data(ii+1,jj,  kk,  n),  &
          Ext_Data(ii,jj+1,kk,  n), Ext_Data(ii+1,jj+1,kk,  n),  &
          Ext_Data(ii,jj,  kk+1,n), Ext_Data(ii+1,jj,  kk+1,n),  &
          Ext_Data(ii,jj+1,kk+1,n), Ext_Data(ii+1,jj+1,kk+1,n) )

          IF ( data_out_high(ind,n)  >   max_mono )            &
               data_out_high(ind,n) = max_mono
        ELSE

          min_mono = MIN (                                       &
          Ext_Data(ii,jj,  kk,  n), Ext_Data(ii+1,jj,  kk,  n),  &
          Ext_Data(ii,jj+1,kk,  n), Ext_Data(ii+1,jj+1,kk,  n),  &
          Ext_Data(ii,jj,  kk+1,n), Ext_Data(ii+1,jj,  kk+1,n),  &
          Ext_Data(ii,jj+1,kk+1,n), Ext_Data(ii+1,jj+1,kk+1,n) )

          IF ( data_out_high(ind,n)  <   min_mono )           &
               data_out_high(ind,n) = min_mono
        END IF

#endif

  END DO
!$OMP END DO NOWAIT
END DO
!$OMP END PARALLEL

! End of routine.
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Mono_Enforce
