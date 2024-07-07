! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This subroutine performs the Singular Value Decomposition (SVD) of matrix A,
! such as A= U' D V, where ' denotes traspose.
! U (V) is the matrix of the left (right) singular vectors of A and D is a
! diagonal matrix whose values are the singular values of A.
!
! Algorithms are based on those described in Golub and Van Loan (1996):
! Matrix Computations, The John Hopkins University press,
! Baltimore and London, 723pp
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Climate Diagnostics
! [OBS!!! If it is at control/misc,
! does it really belong to Climate diagnostics??]

MODULE SVD

! DrHook modules
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! DrHook variables
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SVD'

CONTAINS

SUBROUTINE givens(x_in,y_in,h_out,c_out,s_out)
!-----------------------------------------------------------------------
! This subroutine calculates givens rotations such as for a vector
! x=[x,y] there is a matrix G composed of sin(theta) and cos(theta) such as
!
!        |  c s | | x |   | x^2 + y^2 |
!  G x = | -s c | | y | = |     0     |
!
! Inputs arguments:
!    x,y: Coordinates of the vector
! Output argumens:
!      c_out: cosine of theta
!      s_out: sine of theta
!      h: Norm of input vector

IMPLICIT NONE
! Argument variables
REAL, INTENT(IN) :: x_in ! First coordinate of input vector x
REAL, INTENT(IN) :: y_in ! Second coordinate of input vector x

REAL, INTENT(OUT) :: h_out ! Norm of the input x vector
REAL, INTENT(OUT) :: c_out ! Cosine value of G matrix
REAL, INTENT(OUT) :: s_out ! Sine value of G matrix

! Local variables
REAL  :: absx, absy ! Absolute values of x_in and y_in.

CHARACTER(LEN=*), PARAMETER  :: RoutineName='GIVENS'

REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

absx=ABS(x_in)
absy=ABS(y_in)

! In order to avoid dividing by zero while computing the norm
! the largest value between absy and absx is placed in the denominator.
IF (absy > absx) THEN
  ! compute norm (also known as radius)
  h_out = absy*SQRT(1.0+(absx/absy)**2)
ELSE
  h_out = absx*SQRT(1.0+(absy/absx)**2)
END IF

! Get cosine
  c_out = x_in / h_out
! Get sine
  s_out = y_in / h_out

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE givens

SUBROUTINE householder(dims,x_in,v_out,beta)
!-----------------------------------------------------------------------
! This Subroutine returns the Household vector v_out for a given
! vector  x_in (whose dimensions are dims). output vector v_out satisfies
! P = I - beta v v' such as P x_in = ||x_in|| e_1 where e is the basis vector.

! Input arguments:
!   dims: Dimension of the vector
!   x_in: input vector
! Output arguments:
!   v_out: Householder vector of x_in.
!   beta: Coefficient for the Householder matrix

IMPLICIT NONE
! Argument variables
INTEGER, INTENT(IN) :: dims        ! Dimension of the x vector
REAL,    INTENT(IN) :: x_in(dims)  ! Input vector
REAL,    INTENT(OUT):: beta        ! Coefficient for the Householder
                                   ! matrix P = I -beta v v'
REAL,    INTENT(OUT):: v_out(dims) ! householder vector of x_in

! Local variables
REAL :: sigma ! Norm of input vector (from the second coordinate)
REAL :: mu    ! Norm of householder vector

CHARACTER(LEN=*), PARAMETER  :: RoutineName='HOUSEHOLDER'

REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Compute Norm of input vector (from 2nd coord)
sigma = DOT_PRODUCT( x_in(2:dims),x_in(2:dims) )

! Define the values of v_out from x_in
v_out(1) = 1.0
v_out(2:dims) = x_in(2:dims)
beta = 0.0

! Compute first value of v_out
IF (sigma /= 0.0) THEN
  mu = SQRT( x_in(1)*x_in(1) + sigma)
  IF ( x_in(1) <= 0 ) THEN
    v_out(1) = x_in(1) - mu
  ELSE
    v_out(1) = -sigma/ ( x_in(1) + mu)
  END IF
END IF

! Get beta and normalize householder vector
beta = 2.0* v_out(1)*v_out(1)/( sigma + v_out(1)*v_out(1))
v_out(1:dims)=v_out(1:dims)/v_out(1)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE householder

SUBROUTINE bi_diagonal(m, n, A, ug, vg, d, f)
!-----------------------------------------------------------------------
! This subroutine computes and returns the bi-diagonal matrix of A(m,n)
! as well as the U(m,n), the product of Householder matrices of columns
! and V(n,n), the product of householder matrices of rows.
!
! The algebraic operation is B = U A V,
! where U = U_n U_n-1 ... U_1 and V = V_n-2 ... V_1
!
! For a 5x3 matrix, B has the form:
!
!      | d1  f1    |
!      |     d2 f2 |
!  B = |        d3 |
!      |           |
!      |           |
!
! Input arguments:
!  m,n dimensions of the matrix (No of rows and columns respectively)
!
! Input/Output arguments
!  A: The intial A matrix as input and final diagonalized matrix
!  as output.
!  Ug: Product of householder matrices of the left side
!  Vg: Product of householder matrices of the right side
!  d: A vector of the values in the diagonal of B
!  f: A vector of values in the upper diagonal of B

IMPLICIT NONE
! Argument variables
INTEGER, INTENT(IN) :: m, n ! Dimensions of input matrix A

REAL, INTENT(INOUT) ::                                                  &
      A(m,n)
     ! Input matrix, is transformed into its bi-diagonal matrix

REAL, INTENT(OUT) ::                                                    &
      ug(m,n)                                                           &
      ! Product of Householder matrices on the left side
,     vg(n,n)                                                           &
      ! Product of Householder matrices on the right side
,     d(n)                                                              &
      ! Vector of the values in the diagonal of B
,     f(n)
      ! Vector of the values in the upper diagonal of B

! Local variables
INTEGER ::                                                              &
      i,j,k,l                                                           &
    ! Integers for DO loops
,     n_min
    ! Minimun between m and n

REAL  ::                                                                &
      beta1(n), beta2(n)                                                &
     ! Vector to store the beta values for the Householder matrices
     ! for columns (beta1) and rows (beta2) of A
,     v1(m,n)                                                           &
     ! Matrix to store the Householder vectors of columns of A
,     v2(n,n)                                                           &
     ! Matrix to store the Householder vectors of rows of A
,     p(n,n)                                                            &
     ! Matrix to store the intermediate values of Householder matrix
     ! for columns
,     q(MAX(m,n),MAX(m,n))
     ! Matrix to store the intermediate values of Householder matrix
     ! for vectors

CHARACTER(LEN=*), PARAMETER  :: RoutineName='BI_DIAGONAL'

REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialize matrices
v1    = 0.0
v2    = 0.0
p     = 0.0
q     = 0.0
vg    = 0.0
ug    = 0.0
f     = 0.0
d     = 0.0

! Compute minimum between M and N
n_min = MIN(m,n)

DO j=1,n_min
  ! Do Householder projection for the columns of A such as
  ! A(2)= U(1)A = { ||a1||e1, a2(2),a3(2), ..., a2(n) }
  ! A(3)= U(2) A(2) = { ||a1||e1, y e1 + ||a2||e2 ,a3(3), ..., a2(3) }
  !
  !
  !                                              | x  x ... x x |
  !                                              | 0  x ... x x |
  ! A(n_min) = U(n_min) U(n_min-1) ... U(1) A =  |    .......   |
  !                                              | 0  0     x x |
  !                                              | 0  0       x |
  !
  !where X() are matrices of dimensions m x n and a1 are the
  ! vectors of A(N)

  ! Get Householder vector v1 of A(j:m,j)
  CALL householder(m-j+1,A(j:m,j),v1(j:m,j),beta1(j))

  ! Compute A2 = U A where U is the Householder Matrix of v1(j),
  ! this is done in two steps (UA)_0 = v1 A
  DO l=j,n
    DO i=j,m
      p(l,j)=p(l,j)+v1(i,j)*A(i,l)
    END DO
  END DO
 ! and UA = -beta1 v1 (UA)_0
  DO l=j,n
    DO i=j,m
      A(i,l) = A(i,l) - beta1(j)*v1(i,j)*p(l,j)
    END DO
  END DO

  IF ( j <= n-2 ) THEN
    ! Do Householder projection for the rows of A such as
    ! A= AV(1) ... V(N-2) such as
    !
    !                | x x 0 0 ... 0 0 |
    !                |   x x x ... x x |
    ! A(1) = AV(1) = |     .......     |
    !                |             x x |
    !                |               X !
    !

    ! Get householder vector from A(j,j+1:n)
    CALL householder(n-j,A(j,j+1:n),v2(j+1:n,j),beta2(j))

    ! Compute A2= A V where V is the Householder matrix of v2(j)
    ! Same procedure as befor
    DO l=j,n
      DO i=j,m
        q(i,j) = q(i,j) + A(i,l)*v2(l,j)
      END DO
    END DO

    DO l=j,n
      DO i=j,m
        A(i,l) = A(i,l) - beta2(j)*q(i,j)*v2(l,j)
      END DO
    END DO

  END IF ! End if over left Householder projection

END DO ! end loop over n_min

! Save B's upper diagonal values in f
DO i=1,n-1
  f(i+1) = A(i,i+1)
END DO
! In case n>m also save the n+1 value
IF (n>m) f(n_min+1) = A(n,n+1)

! Save B's diagonal values in d, and compute minimal difference
DO i=1,n
  d(i) = A(i,i)
END DO

! Reconstruct U matrix
p = 0.0
! Build the identity matrix I(n,n)
DO i = 1,n_min
  UG(i,i) = 1.0
END DO

! Loop until n doing U(n) = U(n-1) ... U(2) U(1)
DO j=n_min,1,-1
  ! First do U(n) = v1 U(n-1)
  DO l=j,n
    DO i=j,m
      p(l,j) = p(l,j) + v1(i,j)*ug(i,l)
    END DO
  END DO
  ! Then U(n) = - beta v1 U(n-1)
  DO l=j,n
    DO i=j,m
      ug(i,l) = ug(i,l) - beta1(j)*v1(i,j)*p(l,j)
    END DO
  END DO

END DO

! Reconstruct V matrix
q=0.0
! Build the identity matrix I(n,n)
DO i =1,n
  VG(i,i) = 1.0
END DO

! Loop until n doing V(n) = V(n-1) ... V(2) V(1)
DO j=n,1,-1
  ! First do V(n) = v2 V(n-1)
  DO l=j+1,n
    DO i=j+1,n
      q(l,j)= q(l,j) + v2(i,j)*vg(i,l)
    END DO
  END DO
  ! Then V(n) = - beta v2 V(n-1)
  DO l=j+1,n
    DO i=j+1,n
      vg(i,l) = vg(i,l) - beta2(j)*v2(i,j)*q(l,j)
    END DO
  END DO

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE bi_diagonal

SUBROUTINE svd_decomp(m, n, ug, vg, d, f)
!-----------------------------------------------------------------------
! Calculate the SVD decomposition of a bi-diagonal matrix A using
! the Golub-Reinsch (GR) algorithm.
! A = U D V' where D is a diagonal matrix with the singular values of A.
!
! The GR algorithm splits the input bi-diagonal matrix in three parts,
! such as
!
!
!        | B11         |
!    B = |     B22     |
!        |         B33 |
!
! Where B11(pxp), B22(n-p-q x n-p-q) is upper diagonal and
! B33 (q,q) is diagonal.
! The GR determinies the smallest p and the largest q, and then
! reduces the upper diagonal elements of B22 through a
! Givens rotation (for rows where d = 0.0) or the Golub-Kahan (GK)
! algorithm if d and f /= 0 in a row.
!
! Input variables
!    m,n dimensions of A matrix (rows x columns)
!
! In/Out variables
!    ug: U matrix (already holding the product of householder
!        transformations) is multiplied by the product of
!        givens rotation matrices hat(U).
!    vg: V matrix (already holding the product of householder
!        transformations) is multiplied by the product of
!        givens rotation matrices  hat(V).
!    d: Vector with values of the diagonal of the input matrix
!    f: Vector with the values of the upper diagonal of the input matrix
!       (one row above) e.g:
!
!       | d1 f2    |
!       |    d2 f3 |
!       |       d3 |

IMPLICIT NONE
! Argument variables
INTEGER, INTENT(IN) :: m, n ! dimensions of the input matrix

REAL, INTENT(INOUT) ::                                                  &
      ug(m,n)                                                           &
     ! Unitary matrix with the left singular values
,     vg(n,n)                                                           &
     ! Unitary matrix with the right singular values
,     d(n)                                                              &
     ! Vector of the values in the diagonal of B
,     f(n)
     ! Vector of the values in the upper diagonal of B

! Local variables
INTEGER ::                                                              &
      i,j,k,l                                                           &
     ! indexes for Do loops
,     lstart
     ! Pointer of coordinate where d(k) has f(k) = 0.0

REAL ::                                                                 &
      dk_mu, dk_mu0                                                     &
     ! d(k) - mu, where mu is the eigenvalue of T = B22' x B22,
     ! the dk_mu0 corresponds to the first iteration in the fast
     ! Givens transform
,     x, y, y0, z, h, c, s                                              &
     ! Values of the coordinates and cos(theta), sin(theta) and norm
     ! for Givens rotations
,     c1, c2                                                            &
     ! Save values for u(k,i) and u(k-1,i) when doing UxG, same for V.
,     eps                                                               &
     ! Smallest number for the machine to see 1.0 + eps > 1.0
,     dl                                                                &
     ! Saved value of a coordinate of d vector.
,     anorm
     ! Maximun of the sum of elements in the diagonal
     ! and upper diagonal for same column (d + f )

CHARACTER(LEN=*), PARAMETER  :: RoutineName='SVD_DECOMP'

REAL(KIND=jprb)               :: zhook_handle

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Initialize values
y0=0.0

! Define epsilon to filter f values
eps = EPSILON(1.0)

! Compute the maximun value of d+f
anorm = 0.0
DO i=1,n
  anorm = MAX(anorm,ABS(d(i))+f(i))
END DO

! Do a loop from k=n to 2 so q goes from 0 to n-1.
! The GR runs until q = n (k reaches 1).
DO k=n,2,-1

  ! set to 0.0 values below epsilon
  IF (eps>ABS(d(k))) d(k)=0.0
  IF (eps>ABS(f(k))) f(k)=0.0

  ! Iterate until f(k) is zero
  DO WHILE( ABS(f(k))+anorm > anorm )

    ! Pick up value of d and its coordinate
    dl=d(k)
    lstart=k

    ! Find last element where f is zero. This is p value.
    ! (size of B11 matrix)
    DO l=k,1,-1
      IF (ABS(f(l)) + anorm == anorm) THEN
        dl=d(l)
        lstart=l
      END IF
    END DO

    IF (ABS(d(k-1)) + anorm == anorm) THEN
      ! If d(k-1) is zero. Apply Givens rotations to the k-th
      ! column so f(k) is reduced to zero and B22
      ! is still upper bi-diagonal
      !
      ! Givens rotation is peformed over vector
      ! x = [d(k) f(k)] such as Gx = [h,0]
      ! where G is given matrix and h norm of x.
      ! U need to be updated such as:
      ! UB'= UGB where B'(k,k+1) = 0

      ! Initialize values prior to the Givens rotation
      x    = d(k)
      y    = f(k)

      ! Perform Givens Rotation as
      !    |x|   | h |
      !  G |y| = | 0 |
      CALL givens(x,y,h,c,s)

      ! Then values of the B22 k row become d(k) = h f(k)= 0.0
      d(k) = h
      f(k) = 0.0

      ! Compute UG, only rows k and k-1 are affected
      DO j=1,m
        ! Save u(j,k-1) and ug(j,k) as c1 and c2
        ! before they are re-computed
        c1 = ug(j,k-1)
        c2 = ug(j,k)

        ! Do u(j,k-1) = c ug (j,k-1 ) + s ug(j,k)
        ug(j,k-1) = c1*c + c2*s
        ! Do u(j,k) = -s ug (j,k) + c ug(j,k)
        ug(j,k)   = -(c1*s) + c2*c
      END DO

    ELSE
      ! The GK algorithm is applied to B22 in order to reduce the value
      ! of the elements in the upper diagonal.
      !
      !
      ! The GK algorithm follows the steps:
      ! - as V' A' A V = D^2; we build a triagonal matrix T = A'A, such as
      !
      !        | .............................................. |
      !        | ... d(k-1)^2 + f(k-1)^2       d(k-1) f(k)  ... |
      !    T = | ...      d(k-1) f(k)       d(k)^2 + f(k)^2 ... |
      !        | .............................................. |
      !
      ! - Computes the eigenvalues of the lower-right 2x2 submatrix of T
      ! - For every l column (p<l<k), apply a Givens rotations to
      !   vector x=[d(k) - mu, d(k) f(k)] and update U=UG, and apply
      !   given rotations to x'G' and update V=GV
      !

      ! Compute d(k) - mu, where mu is the eigenvalue of T(2x2)
      dk_mu0=( (d(k-1)-d(k)) * (d(k-1)+d(k)) +                          &
             (f(k-1)-f(k)) * (f(k-1)+f(k)) ) /                          &
             (2.0*f(k)*d(k-1))

      CALL givens(dk_mu0,1.0,h,c,s)

      ! Apply Fast-Givens transformation to get the first dimension of the
      ! vector x such as Gx=[h,0]
      dk_mu =( (dl-d(k))*(dl+d(k)) +                                    &
           f(k) * ( (d(k-1)/(dk_mu0 + SIGN(h,dk_mu0))) - f(k) ) )/dl

      ! Initialize the cos(theta) and sin(theta) for Givens rotation
      c=1.0
      s=1.0

      ! Apply rotations for each column from lstart (or p) to k-1
      DO l=lstart,k-1

        ! Apply Givens rotation to rows,
        ! define elements of the vector.
        z     =  s*f(l+1)
        y     =  c*f(l+1)

        CALL givens(dk_mu,z,h,c,s)

        ! Update elements in the diagonal after Givens rotation
        f(l)  =  h
        dk_mu =  (dl*c) + (y*s)
        y     =   (y*c) - (dl*s)
        z     =  s*d(l+1)
        y0    =  c*d(l+1)

        ! Update V'= GV
        DO i=1,n
          !Save values
          c1 = vg(i,l)
          c2 = vg(i,l+1)
          !apply G matrix product
          vg(i,l)   = (c*c1) + (s*c2)
          vg(i,l+1) = (c*c2) - (s*c1)
        END DO

        ! Perform Givens Rotation for columns
        CALL givens(dk_mu,z,h,c,s)
        d(l) = h
        dk_mu =  (c*y)+(s*y0)
        dl = -(s*y)+(c*y0)

        ! Update U=UG
        DO i=1,m
          !save values
          c1 = ug(i,l)
          c2 = ug(i,l+1)
          !apply G matrix product
          ug(i,l)   = (c*c1) + (s*c2)
          ug(i,l+1) = (c*c2) - (s*c1)
        END DO

      END DO     ! Loop over l (from p to k-1)

      ! reset values of upper-diagonal and d(k) after GK algorithm
      f(lstart) = 0.0
      f(k) = dk_mu
      d(k) = dl

    END IF ! End if for f(k) > 0.0

  END DO   ! End while until f(k) = 0.0
END DO     ! End Loop over k

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE svd_decomp

END MODULE SVD
