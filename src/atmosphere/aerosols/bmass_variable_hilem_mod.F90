! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module including subroutines to inject elevated emissions of biomass
!   burning aerosol when their corresponding injection heights are not
!   the same for all grid cells.
!
! Method:
!   Three subroutines are included in this module:
!   1) calc_interf_z:  Calculate altitude of interfaces for the model
!      theta layers.
!   2) get_injection_fraction:  Get a 3-D array with fraction of bmass emission
!      in each vertical layer 
!      vertical layers for the injection of elevated emissions 
!   3) inject_variable_emiss: Update fresh biomass burning aerosol mmr
!      with elevated emissions within the 2-D varying injection layers.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: aerosols
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------------

MODULE bmass_variable_hilem_mod

USE parkind1,         ONLY: jpim, jprb     ! DrHook
USE yomhook,          ONLY: lhook, dr_hook ! DrHook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='BMASS_VARIABLE_HILEM_MOD'

CONTAINS

!---------------------------------------------------------------------------
! Description:
!   Calculate altitude of interfaces for the model theta layers.
!
! Method:
!   Use information on model dimensions to calculate a 3-dimensional
!   array with the altitude of the interfaces.
!---------------------------------------------------------------------------
SUBROUTINE calc_interf_z (row_length, rows, model_levels, interf_z)

USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN)    :: row_length    ! Model dimensions
INTEGER, INTENT(IN)    :: rows
INTEGER, INTENT(IN)    :: model_levels

! Altitude of interfaces (metres) above ground level
REAL,    INTENT(INOUT) :: interf_z (row_length, rows, 0:model_levels)


! Local variables
INTEGER :: i,j,k   ! loop counter

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER    :: RoutineName='CALC_INTERF_Z'

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k) SHARED(model_levels, rows,      &
!$OMP row_length, interf_z, r_rho_levels, r_theta_levels) 

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels-1
  DO j=1,rows
    DO i=1,row_length
      interf_z(i, j, k) = r_rho_levels(i, j, k+1) - r_theta_levels (i, j, 0)
    END DO
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO j=1,rows
  DO i=1,row_length
    interf_z(i, j, 0) = 0.0
    interf_z(i, j, model_levels) = r_theta_levels(i, j, model_levels) -     &
                                   r_theta_levels (i, j, 0)
  END DO
END DO
!$OMP END DO

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE calc_interf_z


!---------------------------------------------------------------------------
! Description:
!   Get a 3-D array with fraction of bmass emission in each vertical layer.
!
! Method:
!   For each horizontal location identify whether or not there are biomass
!   burning emissions and if they are within a single grid box. 
!   If there are emissions, calculate the total depth over which they are 
!   emitted (i.e. maximum injection height - minimum injection height). 
!   For each box between the maximum and minimum injection height calculate 
!   the fraction of it's depth relative to the total depth.
!   If all the emissions are within a single box, the fraction is 1.0.
!  
!---------------------------------------------------------------------------
SUBROUTINE get_injection_fraction (row_length, rows, model_levels,           &
                                 bmass_hilem_h1, bmass_hilem_h2, interf_z, &
                                 bmass_frac_em)

USE umPrintMgr,             ONLY: umPrint, umMessage
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN)  :: row_length    ! Model dimensions
INTEGER, INTENT(IN)  :: rows
INTEGER, INTENT(IN)  :: model_levels

! Minimum and maximum heights for injection of elevated biomass burning emiss
REAL,    INTENT(IN)  :: bmass_hilem_h1 (row_length,rows)
REAL,    INTENT(IN)  :: bmass_hilem_h2 (row_length,rows)
! Interfaces of theta model layers
REAL,    INTENT(IN)  :: interf_z       (row_length, rows, 0:model_levels)

! 
REAL, INTENT(OUT) :: bmass_frac_em (row_length,rows,model_levels)


! Local variables
INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

INTEGER  :: i, j, k   ! loop counters
INTEGER  :: icode

REAL     :: interf1, interf2, emiss1, emiss2, emission_depth
REAL     :: bmass_frac_total(row_length,rows)

CHARACTER (LEN=errormessagelength) :: cmessage  ! error message
CHARACTER(LEN=*), PARAMETER        :: RoutineName='GET_INJECTION_FRACTION'

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! for every grid box calculate the fraction of the emissions layer
! which lies in that grid box
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,k,interf1,interf2,emiss1,   &
!$OMP emiss2, emission_depth)  DEFAULT(NONE)                               &
!$OMP SHARED(model_levels,rows,row_length,interf_z,bmass_hilem_h1,         &
!$OMP bmass_hilem_h2,bmass_frac_em)
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length

      ! bottom and top heights of this box
      interf1 = interf_z(i,j,k-1)
      interf2 = interf_z(i,j,k)
      
      ! test whether this box is outside the plume - if so
      ! fraction is zero
      IF (interf2 < bmass_hilem_h1(i,j) .OR. interf1 > bmass_hilem_h2(i,j)) THEN
        bmass_frac_em(i,j,k) = 0.0

      ! test if the plume is entirely within this box - in this case 
      ! the fraction for the box is 1. Note that this means that if 
      ! the top and bottom of the plume are zero then the whole
      ! mass will go in level 1 as expected
      ELSE IF (interf1 <= bmass_hilem_h1(i,j) .AND.  &
               interf2 >  bmass_hilem_h2(i,j)) THEN
        bmass_frac_em(i,j,k) = 1.0

      ! if the box does lie partly within the plume, 
      ! calculate the fraction of the plume within the box
      ELSE 

        ! bottom and top heights of the emissions layer
        emiss1 = bmass_hilem_h1(i,j) 
        emiss2 = bmass_hilem_h2(i,j) 
        
        ! to prevent floating overflow make the emission layer at
        ! least 1m deep
       
        IF ((emiss2 - emiss1) < 1.0 ) emiss2 = emiss1 + 1.0

        ! calculate the depth of the layer with 
        ! emisisons
        emission_depth = emiss2 - emiss1
        
        ! now calculate the fraction. If the box lies partly in the 
        ! emissions layer, then take the height of the emission 
        ! layer as the depth
        bmass_frac_em(i,j,k) = ( MIN(interf2, emiss2) -            &
                                 MAX(interf1, emiss1) ) /          &
                                 emission_depth
      
      END IF 

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! check that sum over levels is approximately equal to 1.
bmass_frac_total = SUM(bmass_frac_em, DIM=3)
IF(ANY(ABS(bmass_frac_total(:,:) - 1.0) > 10.0*EPSILON(bmass_frac_total))) THEN 
    icode    = 1
    cmessage = 'Sum over levels of emission layers for ' //            &
                 'biomass burning aerosol emissions not equal to 1'
    CALL ereport (ModuleName//':'//routinename, icode, cmessage)
END IF 

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE get_injection_fraction

!---------------------------------------------------------------------------
! Description:
!   Update fresh biomass burning aerosol mmr with high-level emissions 
!   distributed between the injection layers based on layer depth.
!
! Method:
!   For each horizontal location multiply the high-level emissions
!   by the fraction of emissions to be injected at each vertical layer.
!   Call TRSRCE to update fresh biomass burning aerosol mmr with the 
!   resulting rate.
!---------------------------------------------------------------------------
SUBROUTINE inject_variable_emiss (                                 &
    row_length, rows, model_levels, off_x, off_y, halo_i, halo_j,  &
    timestep, bmass_hilem, bmass_frac_em,                          &
    theta, rho, q, qcl, qcf, exner_rho_levels, bmass_new)

USE trsrce_mod,           ONLY: trsrce

IMPLICIT NONE

! Subroutine arguments

INTEGER, INTENT(IN)  :: row_length                ! model dimensions
INTEGER, INTENT(IN)  :: rows
INTEGER, INTENT(IN)  :: model_levels

INTEGER, INTENT(IN)  :: off_x                     ! size of small halo in i
INTEGER, INTENT(IN)  :: off_y                     ! size of small halo in j
INTEGER, INTENT(IN)  :: halo_i                    ! size of halo in i direction
INTEGER, INTENT(IN)  :: halo_j                    ! size of halo in j direction

REAL, INTENT(IN) :: timestep                      ! atmosphere model timestep

REAL, INTENT(IN) :: bmass_hilem (row_length,rows) ! biomass hi level emissions
! fraction to emit in each layer
REAL, INTENT(IN) :: bmass_frac_em (row_length,rows,model_levels) 

! Additonal arrays needed for the call to TRSRCE
REAL, INTENT(IN) :: theta (1-off_x: row_length+off_x,   1-off_y: rows+off_y,   &
                            model_levels)
REAL, INTENT(IN) :: rho   (1-off_x: row_length+off_x,   1-off_y: rows+off_y,   &
                            model_levels)

REAL, INTENT(IN) :: q     (1-halo_i: row_length+halo_i, 1-halo_j: rows+halo_j, &
                            model_levels)
REAL, INTENT(IN) :: qcl   (1-halo_i: row_length+halo_i, 1-halo_j: rows+halo_j, &
                            model_levels)
REAL, INTENT(IN) :: qcf   (1-halo_i: row_length+halo_i, 1-halo_j: rows+halo_j, &
                            model_levels)

REAL, INTENT(IN) :: exner_rho_levels  (1-off_x: row_length+off_x, &
                                       1-off_y: rows+off_y,       &
                                       model_levels+1)

! mmr fresh biomass burning aerosol
REAL, INTENT(INOUT) :: bmass_new (1-off_x:row_length+off_x, 1-off_y:rows+off_y,&
                                  model_levels)

! Local variables
!
! Min. and max. possible model levels for hilev emiss of biomass burning aerosol
INTEGER :: min_level, max_level 
INTEGER :: i, j, k              ! loop counters

! Quantity of biomass smoke emitted per model level
REAL    ::  biomass_per_level(row_length,rows)


INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INJECT_VARIABLE_EMISS'

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! loop over levels and call trsrce to do actual emissions
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)                &
!$OMP SHARED(model_levels,rows,row_length,biomass_per_level,bmass_hilem,       &
!$OMP bmass_frac_em,off_x,off_y,halo_i,halo_j,theta,q,qcl,qcf,exner_rho_levels,&
!$OMP rho,bmass_new,timestep)
DO k = 1, model_levels
  DO j = 1,rows
    DO i = 1,row_length
      biomass_per_level(i,j) = bmass_hilem(i,j)*bmass_frac_em(i,j,k)
    END DO
  END DO

  ! Inject 2-D highlev emission field within current model level
  CALL trsrce(                                                           &
    rows, row_length, off_x, off_y, halo_i, halo_j,                      &
    theta, q , qcl , qcf , exner_rho_levels, rho,                        &
    bmass_new(:,:,k), biomass_per_level(:,:),                            &
    k, timestep, 1, 1, 0.0 )

END DO
!$OMP END PARALLEL DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN
END SUBROUTINE inject_variable_emiss

END MODULE bmass_variable_hilem_mod
