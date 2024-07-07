! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   A module containing constants/parameters used in the dust scheme
!
MODULE dust_parameters_mod

!
! Description:
!   This module contains declarations for constants and tunable
!   parameters used to diagnose the emission and deposition of
!   mineral dust
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
!
! Code description:
!   Language: Fortran 90
!   This code is written to UM programming standards version 8.2.
!

USE missing_data_mod, ONLY: rmdi, imdi
USE ereport_mod, ONLY: ereport
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!
! Parameters fundamental to the dust scheme
! (Details of dust properties/size divisions etc.)
! =======================================================================

! Prognostic mineral dust aerosol
LOGICAL :: l_dust
! Diagnostic mineral dust aerosol lifting
LOGICAL :: l_dust_diag

! Logical to contain whether using two-bin (true) or six-bin (false) dust
LOGICAL :: l_twobin_dust = .FALSE.
!   (default to six-bin dust)

! Define which dust bins are active
LOGICAL :: l_dust_div1, l_dust_div2, l_dust_div3, l_dust_div4,        &
           l_dust_div5, l_dust_div6

! Define which dust bins are active for LBCs
LOGICAL :: l_dust_div1_lbc, l_dust_div2_lbc, l_dust_div3_lbc, &
           l_dust_div4_lbc, l_dust_div5_lbc, l_dust_div6_lbc

! Number of discrete particle size divisions (bins)
INTEGER :: ndiv              ! number of divisions that can be lifted
                             !  from the surface
INTEGER, PARAMETER :: ndivh=9! number of divisions that can be blown
                             !  horizontally along the surface and contribute
                             !  to the lifting of the 1st NDIV divisions
INTEGER, PARAMETER :: ndivl=6! number of divisions in dust soil ancillaries

INTEGER :: i_dust = imdi                    ! dust scheme setting in namelist
INTEGER, PARAMETER :: i_dust_off = 0        ! No dust
INTEGER, PARAMETER :: i_dust_prognostic = 1 ! Prognostic dust
INTEGER, PARAMETER :: i_dust_diagnostic = 2 ! Diagnostic dust

! The size of particles included in each division
! Note that by using two arrays for the max/min diameter here we can set
! up overlapping divisions. We must take care to make them consistent

REAL, ALLOCATABLE ::  dmax(:)
        ! max diameter of particles in each div.
REAL, ALLOCATABLE ::  dmin(:)
        ! min diameter of particles in each div.
REAL, ALLOCATABLE ::  drep(:)
        ! representative particle diameter
!
! Physical properties of the dust
REAL, PARAMETER :: rhop = 2.65e+3  ! density of a dust particle (quartz)

!
! Parameters used during dust emissions calculations
! =======================================================================
!
! Parameters based on observations/published research
REAL, PARAMETER :: ustd_bas(ndivh) =                       &
       (/ 0.85, 0.72, 0.59, 0.46, 0.33, 0.16, 0.14, 0.18, 0.28 /)
                                      ! impact U*t derived from Bagnold (1941)
REAL, PARAMETER :: horiz_c = 2.61     ! C in horizontal flux calc (White 1979)
REAL, PARAMETER :: vert_a = 13.4      ! A in vertical flux calc(Gillette 1979)
REAL, PARAMETER :: vert_b = -6.0       ! B in vertical flux calc(Gillette 1979)
REAL, PARAMETER :: vert_c = 0.01      ! cgs to si conversion

!
! Input variables needed when using the dust lifting scheme with 1 tile
!
REAL            :: z0_soil            ! Roughness length over bare soil

!
! Tuning parameters - set by namelist
REAL            :: us_am = rmdi       ! ustar correction (multiplic.)
REAL            :: sm_corr = rmdi     ! soil moist. correction factor
REAL            :: horiz_d = rmdi     ! Global tuning param for horizontal
                                      !  (and hence vertical) flux
! Tuning parameters - defined here
REAL, PARAMETER :: us_aa  = 0.0       ! ustar correction (additional)
REAL, PARAMETER :: ust_aa = 0.0       ! ustar_t correction (add)
REAL, PARAMETER :: ust_am = 1.0       ! ustar_t correction (multi.)
!
! Limits used in dust flux calculations
REAL            :: u_s_min = 1.0e-5    ! Minimum val of u_s_std,
                                      !  below which no dust is produced
                                      !  (avoids divide by zero problem)
REAL, PARAMETER :: clay_max = 0.1     ! Max clay fraction.
REAL            :: snowmin = 1.0e-6    ! Min snow depth for dust
REAL            :: h_orog_limit = 150.0! 1/2pk-trough height above which
                                      !  no dust is produced
REAL, PARAMETER :: fland_lim =0.99999 ! No ems if fland<lim as windspeed
                                      !  too high at coastal points
!
! Switch to diagnose vertical flux using a fixed (user definable) size
!  distribution. If set to false, the vertical flux in each bin at a point
!  is poportional to the horizontal flux in that bin at that point.
LOGICAL         :: l_fix_size_dist = .FALSE.

! Namelist components can't be allocatable, so set to six, but in two-bin
! case should only use the first two array positions.
REAL :: size_dist(6) = rmdi

!
! Switch to allow emission of dust on tiles other than the bare soil tile
! 0 allows emission only on bare soil, 1 uses a bare soil radiative frac
! from the LAI, as for surface albedo (scaling by dust_veg_sc_nojules),
! and other methods (2+ could be added later).
INTEGER         :: dust_veg_emiss = imdi
!
! Parameters used during the gravitational settling of dust
! =======================================================================
!

REAL, PARAMETER :: accf = 1.257 ! Cunningham correction factor term A
REAL, PARAMETER :: bccf = 0.4   ! Cunningham correction factor term B
REAL, PARAMETER :: cccf = -1.1  ! Cunningham correction factor term C

!
! The total dust deposition flux (required by ocean biogeochemistry)
! =======================================================================
!
REAL, ALLOCATABLE :: tot_dust_dep_flux(:,:)

!
! Parameters used during the scavenging of dust
! =======================================================================
!
REAL, ALLOCATABLE :: krain_dust(:)
REAL, ALLOCATABLE :: ksnow_dust(:)
!
! Parameters for dust diagnostics
REAL :: pwsdiag_sfc_em = rmdi ! fraction of dust after emission in PWS diags
!
! RUN_Dust namelist via which non-parameter values herein can be set
! =======================================================================
!
NAMELIST /RUN_Dust/                                                      &
     us_am, sm_corr, horiz_d, l_fix_size_dist, dust_veg_emiss,           &
     i_dust, l_twobin_dust, l_dust_div1_lbc, l_dust_div2_lbc,            &
     l_dust_div3_lbc, l_dust_div4_lbc, l_dust_div5_lbc, l_dust_div6_lbc, &
     pwsdiag_sfc_em

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DUST_PARAMETERS_MOD'

CONTAINS
!
! Internal subroutine to check that entries to RUN_dust are consistent
! =======================================================================
!
SUBROUTINE dust_parameters_check( )
!
!   Check that a user-defined emissions size distribution is normalised
!
IMPLICIT NONE

IF (l_fix_size_dist) THEN
  size_dist(:)=size_dist(:)/SUM(size_dist)
END IF

END SUBROUTINE dust_parameters_check


! Internal subroutine to set the default value of size_dist. This may be
! overwritten by setting size_dist in the namelist (hence it's in a separate
! subroutine) or by setting l_twobin_dust. If l_fix_size_dist is set it is
! possible to override size_dist in the namelist as that was the previous
! behaviour, and it's mandatory to use it in the case of l_twobin_dust. If
! l_fix_size_dist is not set the value here is irrelevant as it's not used

SUBROUTINE dust_size_dist_initialise( )

IMPLICIT NONE

size_dist(1:6) =                                                   &
      (/ 0.0005, 0.0049, 0.0299, 0.2329, 0.4839, 0.2479 /)

END SUBROUTINE dust_size_dist_initialise

! Internal subroutine to load the correct parameters into the variables in this
! module depending on whether two- or six- bin dust is being used
SUBROUTINE dust_parameters_load( )

IMPLICIT NONE

INTEGER           :: icode           ! Error code
CHARACTER(LEN=errormessagelength) :: cmessage        ! Error message
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DUST_PARAMETERS_LOAD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Convert i_dust into appropriate logical settings
SELECT CASE(i_dust)
CASE (i_dust_off)
  l_dust      = .FALSE.
  l_dust_diag = .FALSE.
CASE (i_dust_prognostic)
  l_dust      = .TRUE.
  l_dust_diag = .FALSE.
CASE (i_dust_diagnostic)
  l_dust      = .FALSE.
  l_dust_diag = .TRUE.
CASE DEFAULT
  icode = 2
  WRITE(cmessage, '(A,I2,A)') 'i_dust = ', i_dust, ' is not a valid'&
                            // ' setting'
  CALL ereport('DUST_PARAMETERS_LOAD', icode, cmessage)
END SELECT

! Ascertain which dust bins are active
IF (i_dust == i_dust_prognostic) THEN
  l_dust_div1 = .TRUE.
  l_dust_div2 = .TRUE.
  IF (l_twobin_dust) THEN
    l_dust_div3 = .FALSE.
    l_dust_div4 = .FALSE.
    l_dust_div5 = .FALSE.
    l_dust_div6 = .FALSE.
    l_dust_div3_lbc = .FALSE.
    l_dust_div4_lbc = .FALSE.
    l_dust_div5_lbc = .FALSE.
    l_dust_div6_lbc = .FALSE.
  ELSE
    l_dust_div3 = .TRUE.
    l_dust_div4 = .TRUE.
    l_dust_div5 = .TRUE.
    l_dust_div6 = .TRUE.
  END IF
ELSE
  l_dust_div1 = .FALSE.
  l_dust_div2 = .FALSE.
  l_dust_div3 = .FALSE.
  l_dust_div4 = .FALSE.
  l_dust_div5 = .FALSE.
  l_dust_div6 = .FALSE.
  l_dust_div1_lbc = .FALSE.
  l_dust_div2_lbc = .FALSE.
  l_dust_div3_lbc = .FALSE.
  l_dust_div4_lbc = .FALSE.
  l_dust_div5_lbc = .FALSE.
  l_dust_div6_lbc = .FALSE.
END IF

! Mandatory use of fixed size distribution when using two-bin code
IF (l_twobin_dust) l_fix_size_dist = .TRUE.

! Define ndiv
IF (l_twobin_dust) THEN
  ndiv = 2
ELSE
  ndiv = 6
END IF

! Allocate scavenging co-efficients
ALLOCATE(krain_dust(ndiv))
ALLOCATE(ksnow_dust(ndiv))
ALLOCATE(dmax(ndiv))
ALLOCATE(dmin(ndiv))
ALLOCATE(drep(ndiv))


! Set variables depending on how many bins
IF (l_twobin_dust) THEN

  ! scav. coeff. for rain
  krain_dust(1:ndiv) = (/ 4.0e-5, 3.0e-4 /)
  ! scav. coeff. for snow
  ksnow_dust(1:ndiv) = (/ 4.0e-5, 3.0e-4 /)

  dmax(1:ndiv) =  (/ 4.0e-6, 2.0e-5 /)
            ! max diameter of particles in each div.
  dmin(1:ndiv) =  (/ 2.0e-7, 4.0e-6 /)
            ! min diameter of particles in each div.
  drep(1:ndiv) =  (/ 2.0e-6, 8.0e-6 /)
            ! representative particle diameter

  ! Two-bin dust must have this size distribution
  size_dist(1:ndivl) =                                                &
      (/ 0.1800, 0.8200, 0.0000, 0.0000, 0.0000, 0.0000 /)

ELSE ! Six-bin dust

  ! scav. coeff. for rain
  krain_dust(1:ndiv) =                                                &
   (/ 2.0e-5, 2.0e-5, 3.0e-5, 6.0e-5, 4.0e-4, 4.0e-4 /)
  ! scav. coeff. for snow
  ksnow_dust(1:ndiv) =                                                &
   (/ 2.0e-5, 2.0e-5, 3.0e-5, 6.0e-5, 4.0e-4, 4.0e-4 /)


  dmax(1:ndiv) =                                                      &
   (/ 2.0e-7, 6.32456e-7, 2.0e-6, 6.32456e-6, 2.0e-5, 6.32456e-5 /)
            ! max diameter of particles in each div.
  dmin(1:ndiv) =                                                      &
   (/ 6.32456e-8, 2.0e-7, 6.32456e-7, 2.0e-6, 6.32456e-6, 2.0e-5 /)
            ! min diameter of particles in each div.
  drep(1:ndiv) =                                                      &
                 (/ 0.112468e-06, 0.355656e-06, 0.112468e-05,         &
                    0.355656e-05, 0.112468e-04, 0.355656e-04 /)
            ! representative particle diameter


END IF ! l_twobin_dust

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE dust_parameters_load


! internal subroutine to deallocate the dust arrays - provided for consistency,
! not actually called at present
SUBROUTINE dust_parameters_unload( )

IMPLICIT NONE
DEALLOCATE (drep)
DEALLOCATE (dmin)
DEALLOCATE (dmax)
DEALLOCATE (ksnow_dust)
DEALLOCATE (krain_dust)

END SUBROUTINE dust_parameters_unload

SUBROUTINE print_nlist_run_dust()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_DUST'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_dust', &
    src='dust_parameters_mod')

WRITE(lineBuffer,*)' us_am = ',us_am
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*)' sm_corr = ',sm_corr
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*)' horiz_d = ',horiz_d
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*)' l_fix_size_dist = ',l_fix_size_dist
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*)' dust_veg_emiss = ',dust_veg_emiss
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*)' l_twobin_dust = ',l_twobin_dust
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*) ' l_dust_div1_lbc = ',l_dust_div1_lbc
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*) ' l_dust_div2_lbc = ',l_dust_div2_lbc
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*) ' l_dust_div3_lbc = ',l_dust_div3_lbc
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*) ' l_dust_div4_lbc = ',l_dust_div4_lbc
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*) ' l_dust_div5_lbc = ',l_dust_div5_lbc
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,*) ' l_dust_div6_lbc = ',l_dust_div6_lbc
CALL umPrint(lineBuffer,src='dust_parameters_mod')
WRITE(lineBuffer,'(A,F16.4)') ' pwsdiag_sfc_em = ',pwsdiag_sfc_em
CALL umPrint(lineBuffer,src='dust_parameters_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='dust_parameters_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_dust

#if !defined(LFRIC)
SUBROUTINE read_nml_run_dust(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_DUST'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER  :: n_int = 2
INTEGER, PARAMETER  :: n_real = 4
INTEGER, PARAMETER  :: n_log = 8

TYPE my_namelist
  SEQUENCE
  INTEGER :: dust_veg_emiss
  INTEGER :: i_dust
  REAL :: us_am
  REAL :: sm_corr
  REAL :: horiz_d
  REAL :: pwsdiag_sfc_em
  LOGICAL :: l_fix_size_dist
  LOGICAL :: l_twobin_dust
  LOGICAL :: l_dust_div1_lbc
  LOGICAL :: l_dust_div2_lbc
  LOGICAL :: l_dust_div3_lbc
  LOGICAL :: l_dust_div4_lbc
  LOGICAL :: l_dust_div5_lbc
  LOGICAL :: l_dust_div6_lbc
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,      &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Dust, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Dust", iomessage)

  my_nml % dust_veg_emiss  = dust_veg_emiss
  my_nml % i_dust          = i_dust
  my_nml % us_am           = us_am
  my_nml % sm_corr         = sm_corr
  my_nml % horiz_d         = horiz_d
  my_nml % pwsdiag_sfc_em = pwsdiag_sfc_em
  my_nml % l_fix_size_dist = l_fix_size_dist
  my_nml % l_twobin_dust   = l_twobin_dust
  my_nml % l_dust_div1_lbc = l_dust_div1_lbc
  my_nml % l_dust_div2_lbc = l_dust_div2_lbc
  my_nml % l_dust_div3_lbc = l_dust_div3_lbc
  my_nml % l_dust_div4_lbc = l_dust_div4_lbc
  my_nml % l_dust_div5_lbc = l_dust_div5_lbc
  my_nml % l_dust_div6_lbc = l_dust_div6_lbc

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  dust_veg_emiss  = my_nml % dust_veg_emiss
  i_dust          = my_nml % i_dust
  us_am           = my_nml % us_am
  sm_corr         = my_nml % sm_corr
  horiz_d         = my_nml % horiz_d
  pwsdiag_sfc_em = my_nml % pwsdiag_sfc_em
  l_fix_size_dist = my_nml % l_fix_size_dist
  l_twobin_dust   = my_nml % l_twobin_dust
  l_dust_div1_lbc =  my_nml % l_dust_div1_lbc
  l_dust_div2_lbc =  my_nml % l_dust_div2_lbc
  l_dust_div3_lbc =  my_nml % l_dust_div3_lbc
  l_dust_div4_lbc =  my_nml % l_dust_div4_lbc
  l_dust_div5_lbc =  my_nml % l_dust_div5_lbc
  l_dust_div6_lbc =  my_nml % l_dust_div6_lbc

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_dust
#endif

END MODULE dust_parameters_mod

