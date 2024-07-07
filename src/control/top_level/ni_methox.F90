! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Calculate chemical changes to water vapour due to methane oxidation
!   and hydrogen photolysis, using the method used at ECMWF
!   (Untch et al (ECMWF Newsletter No 82, winter 1998/99, pp 2-8) and
!   Simmons (pers. comm.)).
!
! Method: Follows (Untch et al (ECMWF Newsletter No 82, winter
!   1998/99, pp 2-8) and Simmons (pers. comm.)). The model methane
!   mixing ratio is implicit, and derived from the assumption that
!   2 [CH4] + [H2O] = 3.75  ppmm throughout the stratosphere.
!   The methane oxidation and hydrogen photolysis rate coefficients
!   vary only with pressure, and are calculated
!   within the code.
!   For the pressures on each model level, this uses a simple 1d array 
!   assuming a surface pressure of 1000 hPa.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Top Level
!
! Code description:
!   Language: Fortran 95.
!   This code is written to UMDP3 standards

! Subroutine Interface:
SUBROUTINE ni_methox(row_length,rows,                         &
     eta_theta_levels,timestep,stashwork,q_n,q_inc)

USE planet_constants_mod, ONLY: pref
USE conversions_mod, ONLY: pi, rsec_per_day 

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE nlsizes_namelist_mod, ONLY: model_levels

USE level_heights_mod, ONLY: z_top_of_model => z_top_theta

USE science_fixes_mod, ONLY: l_methox_fix

USE um_parvars, ONLY: at_extremity
  

USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                  &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf

USE model_domain_mod, ONLY: model_type, mt_single_column

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
  
!   Model dimensions
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows

!   Height of model levels (in eta coordinates)
REAL, INTENT(IN) :: eta_theta_levels(0:model_levels)

!   Model parameters
REAL, INTENT(IN) :: timestep

!   STASH workspace
REAL, INTENT(INOUT) :: STASHwork(*)

!   Humidity field and its increment
REAL, INTENT(IN) :: q_n(row_length, rows, model_levels)
REAL, INTENT(INOUT) :: q_inc(row_length, rows, model_levels)

! Local variables
INTEGER :: i, j, k ! loop counters

!   Maximum of 2CH4+H2O used to imply methane amount
REAL, PARAMETER :: max2MpW=3.75e-06

REAL :: ak1(model_levels)   ! methane oxidation coefficient
REAL :: ak2(model_levels)   ! hydrogen photolysis coefficient
REAL :: press(model_levels) ! array of assumed pressures

!   Empirical coefficients used in calculation of ak1 and ak2
REAL, PARAMETER :: coeff1_k1=19.0
REAL, PARAMETER :: coeff2_k1=10.0
REAL, PARAMETER :: coeff3_k1=20.0
REAL, PARAMETER :: coeff4_k1=100.0
REAL, PARAMETER :: coeff1_k2=0.3333
REAL, PARAMETER :: coeff2_k2=0.01 
REAL, PARAMETER :: coeff3_k2=3.0
REAL, PARAMETER :: coeff4_k2=100.0
REAL, PARAMETER :: coeff5_k2=0.005
REAL, PARAMETER :: coeff6_k2=0.01

!   Pressure boundaries in ak1 function
REAL, PARAMETER :: p1_k1=50.0    
REAL, PARAMETER :: p2_k1=10000.0 

!   Pressure boundaries in ak2 function
REAL, PARAMETER :: p1_k2=0.1
REAL, PARAMETER :: p2_k2=20.0

!   Scale height (in m) for calculating pressure 
REAL, PARAMETER :: scale_ht=7000.0 

REAL :: alp1 ! used in calculation of ak1
REAL :: alp2 ! used in calculation of ak2
  
REAL :: tau1 ! chemical timescale for methane oxidation
REAL :: tau2 ! chemical timescale for photolysis

! Local variables for stash
REAL, ALLOCATABLE :: q_incr_diagnostic(:, :, :)
INTEGER :: icode, item, im_index
INTEGER, PARAMETER :: sect = 4 ! making use of LSP section

!   Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0 ! DrHook tracing entry
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1 ! DrHook tracing exit
REAL(KIND=jprb)               :: zhook_handle  ! DrHook tracing

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: RoutineName='NI_METHOX'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)


SELECT CASE (model_type)
  CASE DEFAULT
    IF (sf(0,sect)) THEN  ! STASHflag set
      IF (sf(982,sect)) THEN
        ! Save q_inc before updating
        ALLOCATE(q_incr_diagnostic(row_length, rows, model_levels))
        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              q_incr_diagnostic(i,j,k) = q_inc(i,j,k)
            END DO ! i
          END DO ! j
        END DO ! k
      END IF
    END IF

  CASE (mt_single_column)
END SELECT ! model_type

! Calculate ak1 and ak2 coefficients for methane oxidation.
!   Follows the approach of Simmons (pers. comm). Analytical
!   functions are calculated that approximately match chemical
!   lifetime plots shown in Brasseur and Solomon (1986) (the rate
!   coefficient is the reciprocal of the chemical lifetime).
!   It is assumed that the rate coefficients only vary with pressure.
alp1 = coeff1_k1 * LOG(coeff2_k1)/(LOG(coeff3_k1))**4
alp2 = LOG(coeff1_k2+coeff2_k2)

! Calculate simple 1d pressure on model levels
!   Note that prior to vn9.2, this calculation was wrong as the
!   array eta_theta_levels was incorrectly indexed on declaration
!   This incorrect behaviour is currently maintained by default, but
!   is fixed by using the temporary logical l_methox_fix. 
DO j = 1, model_levels

  IF (l_methox_fix) THEN
    press(j) = pref*EXP(-(z_top_of_model*eta_theta_levels(j)/scale_ht))
  ELSE
    press(j) = pref*EXP(-(z_top_of_model*eta_theta_levels(j-1)/scale_ht))
  END IF ! l_methox_fix

  ! Calculate k1 coefficient
  IF (press(j) <= p1_k1) THEN
    ak1(j) = 1/(rsec_per_day*coeff4_k1)
  END IF
  IF (press(j) > p1_k1 .AND. press(j) < p2_k1) THEN
    tau1 = coeff4_k1                                                    &
            * (1+alp1*((LOG(press(j)/p1_k1))**4/LOG(p2_k1/press(j))))
    ak1(j) = 1/(rsec_per_day*tau1)
  END IF
  IF (press(j) >= p2_k1) THEN
    ak1(j)=0.0
  END IF

  ! Calculate k2 coefficient
  IF (press(j) <= p1_k2) THEN
    ak2(j) = 1/(rsec_per_day*coeff3_k2)
  END IF
  IF (press(j) > p1_k2 .AND. press(j) < p2_k2) THEN
    tau2 = 1 / (EXP(alp2-0.5*(LOG(coeff4_k2)+alp2)                      &
                    * (1 + COS(pi*LOG(press(j)/p2_k2)                   &
                        / LOG(coeff5_k2))))-coeff6_k2)
    ak2(j)=1/(rsec_per_day*tau2)
  END IF
  IF (press(j) >= p2_k2) THEN
    ak2(j)=0.0
  END IF
    
END DO

! Add on q increments
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k) SHARED(model_levels,rows,  &
!$OMP row_length,q_n,q_inc,ak1,ak2,timestep)
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      IF (q_n(i,j,k) >  0.0 .AND. q_n(i,j,k) <  max2MpW) THEN
        q_inc(i,j,k) = q_inc(i,j,k) +                                   &
             (ak1(k)*(max2MpW-q_n(i,j,k))-ak2(k)*q_n(i,j,k))*timestep
      END IF
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL
    

! ----------------------------------------------------------------------
! Write out STASH diagnostics
! ----------------------------------------------------------------------

SELECT CASE (model_type)
  CASE DEFAULT

    IF (sf(0,sect)) THEN ! diagnostics requested this timestep

      icode = 0 ! Initialise error status
      im_index = 1

      item = 982
      IF (icode <= 0 .AND. sf(item,sect)) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              q_incr_diagnostic(i,j,k) = q_inc(i,j,k) - q_incr_diagnostic(i,j,k)
            END DO !i
          END DO  !j
        END DO   !k

        ! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
             q_incr_diagnostic,                                        &
             row_length,rows,model_levels,0,0,0,0, at_extremity,       &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
             stash_levels,num_stash_levels+1,                          &
             atmos_im,sect,item,                                       &
             icode,cmessage)
        DEALLOCATE(q_incr_diagnostic)

        IF (icode >  0) THEN
          cmessage=": error in copydiag_3d(item 982)"//cmessage
        END IF
      END IF

      IF (icode /= 0) THEN
        CALL Ereport(RoutineName,icode,Cmessage)
      END IF

    END IF  ! on sf(0,sect)

  CASE (mt_single_column)
END SELECT ! model_type

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ni_methox
