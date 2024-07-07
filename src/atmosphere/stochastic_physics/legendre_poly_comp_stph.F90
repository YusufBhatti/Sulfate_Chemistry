! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This routine computes the Associate Legendre Polynomials and populates
! a matrix Ymn with their values for different latitudes.
! 
! It employs recurrence relationships from Y(0,0)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics

MODULE legendre_poly_comp_stph_mod

IMPLICIT NONE

LOGICAL, PARAMETER :: create_legendre = .TRUE.
LOGICAL, PARAMETER :: destroy_legendre = .FALSE.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LEGENDRE_POLY_COMP_STPH_MOD'

CONTAINS

SUBROUTINE legendre_poly_comp_stph(global_rows, delta_phi, l_create)

 USE stochastic_physics_run_mod, ONLY: l_x_eq_sin_x, stph_n2, Ymn
 
 USE conversions_mod, ONLY: pi

 USE yomhook, ONLY: lhook, dr_hook
 USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

REAL, INTENT(IN) ::                                                     &
       delta_phi                                                             
      ! grid latitude  spacing in radians
INTEGER, INTENT(IN) ::                                                  &
       global_rows
      ! No of latitude points

LOGICAL, INTENT(IN) ::                                                  &
       l_create
      ! Flag to create or deallocate array

INTEGER ::                                                              &
       ntop                                                            
      ! Truncation (stph_n1~lower and stph_n2~upper) limit.

REAL ::                                                                 &
    sin_v(global_rows)                                                  &
    ! array for sine to build up Fm matrix             
,   cos_v(global_rows) 

REAL ::                                                                 &
    SH_fac                                                              &
    ! Factor for SH before j (lat loop) 
,   SH_lm1                                                                
    ! Factor for the Y(l-1,m) term


INTEGER :: l,m,j                                                        &
    ! Indexes for DO loops
,  lp1, lm1
    ! l+1 and l-1 for indexing arrays



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEGENDRE_POLY_COMP_STPH'

! ++++++++++++++++++++++++++++++++++++++++++++++++++++
! End of variable declaration start of the subroutine
! +++++++++++++++++++++++++++++++++++++++++++++++++++

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! +++++++++++++++++++++++++++++++++++++++++++++++++++
!                      Create Ymn 
IF (l_create) THEN 

  ! If linear approximation requested
  IF (l_x_eq_sin_x) THEN 
    ! Create cos_v using the linear approximation cos(x) = 2*x/Nlat-1
    ! where x=0,..,Nlat
    DO j=1,global_rows
      cos_v(j)= 2.*REAL(j)/REAL(global_rows) - 1.0
    END DO  
 
    ! Create sin_v using the value of cos_v
    DO j=1,global_rows
      sin_v(j)=SQRT( 1.0 - cos_v(j)**2) * (-1.0)
    END DO  
  
  ELSE 

    ! Put (-1) for consistency with the former method.
    ! Remove it for future version that don't require a 
    ! field to field comparison of the FP.

    ! Create sin_v: Use Cos (lat-Pi/2) to remap the latitude
    DO j=1,global_rows
      sin_v(j)=SIN(delta_phi*(j-0.5))*(-1.0)
    END DO

    ! Create sin_v -> mind the '-' in the sin(phi -Pi/2) = - cos(phi)
    DO j=1,global_rows
      cos_v(j)=COS( delta_phi*(j-0.5))*(-1.0)
    END DO
    
  END IF ! end if over creation of sin and cos

  ! ----------------- Declaration of Ymn matrix ------------------------
  ! Set maximun wavenumber
  ntop = stph_n2

  IF (.NOT. ALLOCATED(Ymn)) THEN
    ALLOCATE(Ymn( global_rows,0:ntop,0:ntop) )
  END IF

  ! ---------------- COMPUTATION OF LEGENDRE ASSOCIATE FUNCTIONS ----
  SH_fac = 0.5 / SQRT(pi)
  DO j=1,global_rows
    Ymn(j,0,0) = SH_fac
  END DO

  ! Do Y(l+1,l+1) from Y(l,l)
  DO l=0,ntop-1
    lp1=l+1

    SH_fac = (-1.)*SQRT( REAL(2*l + 3)/REAL(2*l+2) )
    DO j=1,global_rows
      Ymn(j,lp1,lp1) = SH_fac * sin_v(j) * Ymn(j,l,l)
    END DO
  END DO  

  ! Do Y(l+1,l) from Y(l,l)
  DO l=0,ntop-1 
    lp1=l+1

    SH_fac = SQRT( REAL(2*l + 3))
    DO j=1,global_rows
      Ymn(j,lp1,l) = SH_fac * cos_v(j) * Ymn(j,l,l)
    END DO
  END DO  

  ! Do Y(l+1,0) from Y(l,0) and Y(l-1,0)
  DO l=1,ntop-1 
    lp1=l+1
    lm1=l-1
    SH_fac = 1.0 / REAL(l+1) * SQRT( REAL((2*l + 3)*(2*l+1)) )

    SH_lm1 = (-1.0)*REAL(l) / REAL(l+1) *                               &
                    SQRT( REAL(2*l + 3) / REAL(2*l - 1) )

    DO j=1,global_rows
      Ymn(j,lp1,0) = SH_fac * cos_v(j)*Ymn(j,l,0) +                     &
                     SH_lm1 * Ymn(j,lm1,0)
    END DO   

  END DO  

  ! Do Y(l+1,m) from Y(l,m) and Y(l-1,m)
  DO m=1,ntop-2
    ! Do it from m+1 to ntop (but mind the l+1 of the lhs)
    DO l=m+1,ntop-1
      lp1= l + 1
      lm1= l - 1

      SH_fac = SQRT( (4.* REAL((l+1)*(l+1)) -1.)  /                     &   
                      (REAL((l+1)*(l+1)) - REAL(m*m) ) )
                       

      SH_lm1 = SH_fac *                                                 &
               SQRT( (REAL(l*l) - REAL(m*m)) / (4.* REAL(l*l) -1.) )
   
      DO j=1,global_rows
        Ymn(j,lp1,m) = SH_fac* cos_v(j) * Ymn(j,l,m) -                  &
                                SH_lm1 * Ymn(j,lm1,m) 
      END DO   
 
    END DO
  END DO  

  ! Set Ymn to zero for m > n
  DO m=1,ntop
    DO l=0,m-1
      DO j=1,global_rows
        Ymn(j,l,m) = 0.0
      END DO
    END DO  
  END DO

ELSE 
! +++++++++++++++++++++++++++++++++++++++++++++++++++
!                      Deallocate Ymn 
  IF (ALLOCATED(Ymn)) DEALLOCATE(Ymn)
  
END IF  

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE legendre_poly_comp_stph

END MODULE legendre_poly_comp_stph_mod
