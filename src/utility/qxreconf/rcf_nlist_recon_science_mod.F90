! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_nlist_recon_science_mod

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

USE missing_data_mod, ONLY: imdi, rmdi

! Description:
!   Defines variables for the RECON_SCIENCE namelist.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


IMPLICIT NONE

PRIVATE

! make the subroutines public
PUBLIC :: read_nml_recon_science,                                            &
          print_nlist_recon_science,                                         &
          check_nml_recon_science

! make Namelist elements public
PUBLIC :: w_zero_start,                w_zero_end,                           &
          q_min,                       snow_landice_min,                     &
          snow_icefree_max,            coast_adj_method,                     &
          use_smc_stress,              polar_check,                          &
          l_adj_t_soil,                l_use_zero_frac_tile_temp,            &
          l_canopy_snow_throughfall,   l_snow_tile_gbm_ancil_fix,            &
          l_force_relayer,             l_regularize_landice,                 &
          l_regularize_landice_all_soil_pts,                                 &
          l_regularize_landice_all_lice_pts

! make useful parameters public
PUBLIC :: coast_adj_standard_method,                                         &
          coast_adj_spiral_method,                                           &
          coast_adj_circle_method

CHARACTER(LEN=*), PARAMETER :: ModuleName='RCF_NLIST_RECON_SCIENCE_MOD'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

!----------------------------------------------------
! Items set from RECON_SCIENCE namelist
!----------------------------------------------------

! Reset w-components of wind to zero - start level
INTEGER       :: w_zero_start = imdi

! Reset w-components of wind to zero - end level
INTEGER       :: w_zero_end = imdi

! Set the minimum value for specific humidity
! For ENDGame prognostics, this is the minimum vapour mixing ratio
REAL          :: q_min = rmdi

! Method for unresolved points, previously lspiral_s (T=>2, F=>1, new method=>3)
INTEGER       :: coast_adj_method = imdi

! Value of coast_adj_method if standard method is to be used
INTEGER, PARAMETER :: coast_adj_standard_method = 1

! Value of coast_adj_method if spiral method is to be used
INTEGER, PARAMETER :: coast_adj_spiral_method   = 2

! Value of coast_adj_method if spiral circle method is to be used
INTEGER, PARAMETER :: coast_adj_circle_method   = 3

! Use soil moisture stress for interpolating soil moisture
LOGICAL       :: use_smc_stress = .FALSE.

! Check output dump for non-uniform polar rows
LOGICAL       :: polar_check = .FALSE.

! Perform canopy snow throughfall
LOGICAL       :: l_canopy_snow_throughfall = .FALSE.

! Force relayering of the snowpack to deal with changes to
! the structure of the snowpack without a change in the
! number of layers
LOGICAL       :: l_force_relayer = .FALSE.

! Regularize snow fields if land ice is reset as ice-free
! and vice versa
LOGICAL       :: l_regularize_landice                 = .FALSE.
! Apply to snow_icefree_max to all soil points
LOGICAL       :: l_regularize_landice_all_soil_pts    = .FALSE.
! Apply to snow_icefree_max to all soil points
LOGICAL       :: l_regularize_landice_all_lice_pts    = .FALSE.
! Minimum permitted snow mass on land ice (kgm-2)
REAL          :: snow_landice_min     = rmdi

! Maximum permitted snow mass on ice-free tiles (kgm-2)
REAL          :: snow_icefree_max     = rmdi

! Adjust soil temperature for newly added hills and valleys
LOGICAL       :: l_adj_t_soil = .FALSE.

! Use value of tstar from tiles with zero fraction in input dump
! for interpolation and output dump
LOGICAL       :: l_use_zero_frac_tile_temp = .FALSE.

! Switch to prevent snow analysis being used until trials are assessed.
! Triggers snow on tile to be set from GBM snow when the latter is read from
! an ancillary.
LOGICAL       :: l_snow_tile_gbm_ancil_fix = .FALSE.

!----------------------------------------------------
! Namelist RECON_SCIENCE.
!----------------------------------------------------

NAMELIST /recon_science/ q_min, w_zero_start, w_zero_end,  &
 use_smc_stress, polar_check, l_canopy_snow_throughfall,   &
 l_force_relayer, coast_adj_method,                        &
 l_regularize_landice, snow_landice_min, snow_icefree_max, &
 l_regularize_landice_all_lice_pts,                        &
 l_regularize_landice_all_soil_pts,                        &
 l_adj_t_soil, l_use_zero_frac_tile_temp,                  &
 l_snow_tile_gbm_ancil_fix


CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE print_nlist_recon_science()
USE umPrintMgr, ONLY: umPrint

IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
CHARACTER (LEN=*), PARAMETER  :: RoutineName='PRINT_NLIST_RECON_SCIENCE'

! Dr Hook
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL umPrint('Contents of namelist recon_science',         &
    src='rcf_nlist_recon_science_mod')

WRITE(lineBuffer,'(A,I0)')' coast_adj_method = ',coast_adj_method
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,E12.5)')' Q_Min = ',Q_Min
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,I0)')' W_Zero_Start = ',W_Zero_Start
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,I0)')' W_Zero_End = ',W_Zero_End
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,L1)')' Use_SMC_Stress = ',Use_SMC_Stress
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,L1)')' polar_check = ',polar_check
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,L1)')' l_canopy_snow_throughfall = ',   &
   l_canopy_snow_throughfall
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,L1)')' l_force_relayer = ',           &
   l_force_relayer
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,L1)')' l_regularize_landice = ',      &
   l_regularize_landice
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,L1)')' l_regularize_landice_all_soil_pts = ',   &
   l_regularize_landice_all_soil_pts
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,L1)')' l_regularize_landice_all_lice_pts = ',   &
   l_regularize_landice_all_lice_pts
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,E12.5)')' snow_landice_min = ',       &
   snow_landice_min
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,E12.5)')' snow_icefree_max = ',       &
   snow_icefree_max
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,L1)')' l_adj_t_soil = ',l_adj_t_soil
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,L1)')' l_use_zero_frac_tile_temp = ',   &
     l_use_zero_frac_tile_temp
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')
WRITE(lineBuffer,'(A,L1)')' l_snow_tile_gbm_ancil_fix = ', &
   l_snow_tile_gbm_ancil_fix
CALL umPrint(lineBuffer,src='rcf_nlist_recon_science_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -',    &
    src='rcf_nlist_recon_science_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE print_nlist_recon_science

!-------------------------------------------------------------------------------

SUBROUTINE read_nml_recon_science(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat
USE errormessagelength_mod, ONLY: errormessagelength
USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 3
INTEGER, PARAMETER :: n_real = 3
INTEGER, PARAMETER :: n_log = 10

TYPE my_namelist
  SEQUENCE
  INTEGER :: w_zero_start
  INTEGER :: w_zero_end
  INTEGER :: coast_adj_method
  REAL    :: q_min
  REAL    :: snow_landice_min
  REAL    :: snow_icefree_max
  LOGICAL :: use_smc_stress
  LOGICAL :: polar_check
  LOGICAL :: l_canopy_snow_throughfall
  LOGICAL :: l_adj_t_soil
  LOGICAL :: l_use_zero_frac_tile_temp
  LOGICAL :: l_snow_tile_gbm_ancil_fix
  LOGICAL :: l_regularize_landice
  LOGICAL :: l_regularize_landice_all_soil_pts
  LOGICAL :: l_regularize_landice_all_lice_pts
  LOGICAL :: l_force_relayer

END TYPE my_namelist

TYPE (my_namelist) :: my_nml
! Dr Hook
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RECON_SCIENCE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log )

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=recon_science, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist recon_science", iomessage)

  my_nml % w_zero_start              = w_zero_start
  my_nml % w_zero_end                = w_zero_end
  my_nml % coast_adj_method          = coast_adj_method
  my_nml % q_min                     = q_min
  my_nml % use_smc_stress            = use_smc_stress
  my_nml % polar_check               = polar_check
  my_nml % l_canopy_snow_throughfall = l_canopy_snow_throughfall
  my_nml % l_force_relayer           = l_force_relayer
  my_nml % l_regularize_landice      = l_regularize_landice
  my_nml % l_regularize_landice_all_soil_pts = l_regularize_landice_all_soil_pts
  my_nml % l_regularize_landice_all_lice_pts = l_regularize_landice_all_lice_pts
  my_nml % snow_landice_min          = snow_landice_min
  my_nml % snow_icefree_max          = snow_icefree_max
  my_nml % l_adj_t_soil              = l_adj_t_soil
  my_nml % l_use_zero_frac_tile_temp = l_use_zero_frac_tile_temp
  my_nml % l_snow_tile_gbm_ancil_fix = l_snow_tile_gbm_ancil_fix

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  w_zero_start              = my_nml % w_zero_start
  w_zero_end                = my_nml % w_zero_end
  coast_adj_method          = my_nml % coast_adj_method
  q_min                     = my_nml % q_min
  use_smc_stress            = my_nml % use_smc_stress
  polar_check               = my_nml % polar_check
  l_canopy_snow_throughfall = my_nml % l_canopy_snow_throughfall
  l_force_relayer           = my_nml % l_force_relayer
  l_regularize_landice      = my_nml % l_regularize_landice
  l_regularize_landice_all_soil_pts = my_nml % l_regularize_landice_all_soil_pts
  l_regularize_landice_all_lice_pts = my_nml % l_regularize_landice_all_lice_pts
  snow_landice_min          = my_nml % snow_landice_min
  snow_icefree_max          = my_nml % snow_icefree_max
  l_adj_t_soil              = my_nml % l_adj_t_soil
  l_use_zero_frac_tile_temp = my_nml % l_use_zero_frac_tile_temp
  l_snow_tile_gbm_ancil_fix = my_nml % l_snow_tile_gbm_ancil_fix

END IF

CALL mpl_type_free(mpl_nml_type,icode)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_recon_science

!-------------------------------------------------------------------------------

SUBROUTINE check_nml_recon_science()

USE chk_opts_mod,           ONLY: chk_var, def_src
USE nlsizes_namelist_mod,   ONLY: model_levels
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE jules_snow_mod,         ONLY: nsmax

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER   :: RoutineName='CHECK_NML_RECON_SCIENCE'
! Dr Hook
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=15)             :: range_string
CHARACTER(LEN=5)              :: str_model_levels
CHARACTER(LEN=5)              :: str_w_zero_start
CHARACTER(LEN=5)              :: str_w_zero_end
REAL                          :: rmdi_tol

CHARACTER (LEN=errormessagelength)  :: cmessage
INTEGER                             :: errorstatus

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

WRITE(str_model_levels,'(I0)')  model_levels
WRITE(str_w_zero_start,'(I0)')  w_zero_start
WRITE(str_w_zero_end,  '(I0)')  w_zero_end
range_string = '[1:' // TRIM(str_model_levels) // ']'
rmdi_tol  = TINY(rmdi)

CALL chk_var( w_zero_start, 'w_zero_start', TRIM(range_string),              &
              cmessage='Value must be -1 or between 1 and ' //               &
                        TRIM(str_model_levels) )
CALL chk_var( w_zero_end, 'w_zero_end', TRIM(range_string),                  &
              cmessage='Value must be -1 or between 1 and ' //               &
                        TRIM(str_model_levels) )
CALL chk_var( w_zero_end, 'w_zero_end',                                      &
              '['//str_w_zero_end//'>='//TRIM(str_w_zero_start)//']',        &
              cmessage='must be >= ' // TRIM(str_w_zero_start) )
CALL chk_var( coast_adj_method, 'coast_adj_method',                          &
              [coast_adj_standard_method, coast_adj_spiral_method,           &
              coast_adj_circle_method],                                      &
              cmessage='Value chosen not one of standard_method,'//          &
                      ' spiral_method or circle_method.')
CALL chk_var( q_min, 'q_min', '[0.0:1.0]',  &
              cmessage='q_min must be between 0 and 1.')

IF ( l_regularize_landice ) THEN
  CALL chk_var( snow_landice_min, 'snow_landice_min', '[>=0.0]',             &
                cmessage='Must be greater than, or equal to, 0')
  CALL chk_var( snow_icefree_max, 'snow_icefree_max', '[>=0.0]',             &
                cmessage='Must be greater than, or equal to, 0')
END IF

IF ( l_regularize_landice_all_lice_pts .OR.                                  &
     l_regularize_landice_all_soil_pts ) THEN
  IF ( .NOT. l_regularize_landice .OR. nsmax > 0 ) THEN
    errorstatus = 10
    WRITE(cmessage,'(A,L1,I6)') 'l_regularize_landice_all_lice_pts ' //      &
       'or l_regularize_landice_all_soil_pts only applies when ' //          &
       'l_regularize_landice with nsmax = 0. Currently these are: ',         &
       l_regularize_landice, nsmax
    CALL ereport( routinename, errorstatus, cmessage )
  END IF
END IF

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nml_recon_science

!-------------------------------------------------------------------------------

END MODULE rcf_nlist_recon_science_mod
