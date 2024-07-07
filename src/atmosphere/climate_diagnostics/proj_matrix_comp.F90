! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This subroutine computes Fm matrix for each m in the third coordinate,
! uses Swarztrauber and Spotz 2003
! where Fm= U x U' where U is the projection matrix of the
! singular value decomposition of the Legendre Polynomial matrix
! Pm=U D V' where D is the diagonal matrix with the enigenvalues.
!
! If Hoskins filter is applied then it needs to get D and V matrices
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Climate Diagnostics

MODULE proj_matrix_comp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PROJ_MATRIX_COMP_MOD'

CONTAINS

SUBROUTINE proj_matrix_comp(global_rows, n1, n2, delta_phi, Sm,           &
                            Hos_filter, Fmn)
! Computes Fm matrix for each m in the third coordinate,
! uses Swarztrauber and Spotz 2003 where Fm= U x U' where U is
! the projection matrix of the singular value decomposition of
! the Legendre Polynomial matrix Pm=U D V' where D is the
! diagonal matrix with the enigenvalues.
!
! If Hoskins filter is applied then we need to get D and V matrices

USE svd,             ONLY: bi_diagonal, svd_decomp ! SVN routines
USE track_mod,       ONLY: Ymn                     ! Legendre Polynomials
USE yomhook,         ONLY: lhook, dr_hook          ! DrHook arg.
USE parkind1,        ONLY: jprb, jpim              ! DrHook arg.

IMPLICIT NONE
! Argument variables
INTEGER, INTENT(IN) :: global_rows   ! No of latitude points.
INTEGER, INTENT(IN) :: n1            ! Lower truncation wavenumber
INTEGER, INTENT(IN) :: n2            ! Upper truncation wavenumber
REAL,    INTENT(IN) :: delta_phi     ! grid latitude  spacing in radians
REAL,    INTENT(IN) :: Sm            ! Amplitude of the smallest wavenumber
LOGICAL, INTENT(IN) :: Hos_filter    ! Logical to activate Hoskins filter
REAL, INTENT(INOUT) :: Fmn(global_rows+1,global_rows+1,0:n2)

                                     ! for the Hoskins filter
!Local variables
INTEGER :: off_m          ! Maximum between m and n1 to to avoid
                          ! points where m > n
INTEGER :: j,n,m,i        ! indexing for DO loops
REAL ::    K_hos          ! Coefficient for Hoskins filter
REAL ::    S_hos(n2-n1+1) ! All Hoskins filter values for n1:n2

! Variables for each m, changes dimensions so need to be allocatable
REAL , ALLOCATABLE ::                                                     &
       Pm (:,:)                                                           &
          ! Matrix with Legendre Polynomials to perform SH trans
,      u_svd(:,:), v_svd(:,:),d_svd(:),f_svd(:)                           &
          ! Matrices of the SVD of Pm
,      Hos_matrix(:,:)                                                    &
          ! Hoskins filter matrix
,      H_Uprim(:,:)
          ! Matrix with the produc H x U' for optimization purposes

CHARACTER(LEN=*), PARAMETER  :: RoutineName='PROJ_MATRIX_COMP'


! DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Create S_hos array if Hoskins filter is requested
IF (Hos_filter) THEN
  K_hos=LOG(1.0/sm)/(REAL(n2)*(REAL(n2)+1.0))**2

  DO n=n1,n2
    S_hos(n-n1+1)=EXP(K_hos*(REAL(n)*REAL(n+1))**2*(-1.0))
  END DO

END IF

! Loop over m to populate Fmn(m)
DO m=0,n2

  ! Extract m>m columns from the Pmn
  off_m=MAX(m,n1)

  ! Allocate Pm (N-m) x (global_rows +1 )
  IF (.NOT. ALLOCATED(Pm))    ALLOCATE(Pm(global_rows+1,n2-off_m+1))
  ! Allocate the rest of matrices for svd
  IF (.NOT. ALLOCATED(u_svd)) ALLOCATE(u_svd(global_rows+1,n2-off_m+1))
  IF (.NOT. ALLOCATED(v_svd)) ALLOCATE(v_svd(n2-off_m+1,n2-off_m+1))
  IF (.NOT. ALLOCATED(f_svd)) ALLOCATE(f_svd(n2-off_m+1))
  IF (.NOT. ALLOCATED(d_svd)) ALLOCATE(d_svd(n2-off_m+1))

  ! Populate Pm value from the Ymn Associated Legendre
  ! Polynomials from Ymn, stored in the track_mod module
  DO n=off_m,n2
    ! Fill it up with the different latitudes
    DO j=1,global_rows+1
      Pm(j,n-off_m+1)= Ymn(j,n,m)
    END DO
  END DO

  ! Do SVD decomposition, first bi-diagonalize matrix
  CALL bi_diagonal(global_rows+1,n2-off_m+1,Pm,u_svd,v_svd,d_svd,f_svd)
  ! Then decompose using theGolub-Reinsch algorithm
  CALL svd_decomp(global_rows+1,n2-off_m+1,u_svd,v_svd,d_svd,f_svd)

  ! Initialize the Fmn matrix
  DO j=1,global_rows+1
    DO i=1,global_rows+1
      Fmn(i,j,m)=0.0
    END DO
  END DO

  ! Apply Hoskins filter if chosen
  IF (Hos_filter) THEN
    ! Allocate matrices for Hoskins filter
    IF (.NOT. ALLOCATED(Hos_matrix))                                      &
       ALLOCATE(Hos_matrix(n2-off_m+1,n2-off_m+1))
    ! Allocate H_uprim for code optimization
    IF (.NOT. ALLOCATED(H_Uprim))                                         &
       ALLOCATE(H_Uprim(global_rows+1,n2-off_m+1))

    ! Create Hoskins matrix
    DO j=1,n2-off_m+1
      DO i=1,n2-off_m+1
        ! Initialize value at i,j
        Hos_matrix(i,j)= 0.0

        !Populate matrix from the SVD values
        DO n=1,n2-off_m+1
          Hos_matrix(i,j) = Hos_matrix(i,j) +                             &
                            S_hos(n+off_m-n1)*v_svd(n,i)*v_svd(n,j)*      &
                            d_svd(i)/d_svd(j)
        END DO
      END DO
    END DO

    !Initiate H_Uprim matrix
    DO i=1,n2-off_m+1
      ! Fill it up with the different latitudes
      DO j=1,global_rows+1
        H_Uprim(j,i)= 0.0
      END DO
    END DO

    ! Create H x U' matrix ( transpose-> To optimize code)
    DO n=1,n2-off_m+1
      DO i=1,n2-off_m+1
        DO j=1,global_rows+1
          H_Uprim(j,i)= H_Uprim(j,i) + Hos_matrix(i,n) * u_svd(j,n)
        END DO
      END DO
    END DO

    ! Compute Fmn (i,j,m) doing the U x (HxU').
    DO n=1,n2-off_m+1
      DO j=1,global_rows+1
        DO i=1,global_rows+1
          Fmn(i,j,m)= Fmn(i,j,m) +  u_svd(i,n)*H_Uprim(j,n)
        END DO
      END DO
    END DO

    ! Deallocate Hos_matrix and H_Uprim
    IF (ALLOCATED(H_Uprim))    DEALLOCATE(H_Uprim)
    IF (ALLOCATED(Hos_matrix)) DEALLOCATE(Hos_matrix)

  ELSE
    ! If Hos filter is not active then Hos_matrix is an identity matrix
    ! So Fm = U x U'
    DO n=1,n2-off_m+1
      DO j=1,global_rows+1
        DO i=1,global_rows+1
          Fmn(i,j,m)= Fmn(i,j,m) +  u_svd(i,n)*u_svd(j,n)
        END DO
      END DO
    END DO

  END IF ! End loop over Hoskins filter

  ! Dealocate arrays for each m
  IF (ALLOCATED(d_svd)) DEALLOCATE(d_svd)
  IF (ALLOCATED(f_svd)) DEALLOCATE(f_svd)
  IF (ALLOCATED(v_svd)) DEALLOCATE(v_svd)
  IF (ALLOCATED(u_svd)) DEALLOCATE(u_svd)
  IF (ALLOCATED(Pm))    DEALLOCATE(Pm)

END DO ! End loop over m

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE proj_matrix_comp

END MODULE proj_matrix_comp_mod
