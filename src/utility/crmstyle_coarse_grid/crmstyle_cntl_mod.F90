! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Module controlling program crmstyle_coarse_grid

MODULE crmstyle_cntl_mod

USE missing_data_mod, ONLY: rmdi, imdi

IMPLICIT NONE
SAVE

! Description:
!   Module containing logical switches to control the program
!   crmstyle_coarse_grid
!
! Method:
!
!------------------------------------------------------------------------------
!   x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|H x x|A x x|A x x|A x x|A x x|H x x|A x x|A x x|A x x|A x x|x x x x
!      -------------------------------------------------------------
!   x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|x x x x
!      -------------------------------------------------------------
!   x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|x x x x
!      -------------------------------------------------------------
!   x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|x x x x
!      -------------------------------------------------------------
!   x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|A x x|x x x x
!      -------------------------------------------------------------
!   x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x|x x x x
!      |     |     |     |     |     |     |     |     |     |     |
!   x x|S x x|A x x|A x x|A x x|A x x|H x x|A x x|A x x|A x x|A x x|x x x x
!      -------------------------------------------------------------
!   x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
!
! shows bottom of grid
!   x - theta grid points original model grid
!   S - start of area to be processed  - 3rd column, 2 row in.
!   H - Starting positions of coarsest grid (multiple of A grid, 15x15 points)
!   A - starting poistions of other coarser grid. (3x3 points)
!------------------------------------------------------------------------------

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: crmstyle_coarse_grid
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! Declarations:

INTEGER ::                   &
    nx_start  = imdi         &  ! starting gridpoint in EW direction
   ,ny_start  = imdi         &  ! starting gridpoint in NS direction
   ,nres      = imdi         &  ! Number of coarse resolutions required
                                ! max 10 values?
   ,new_res(10)              &  ! New resolutions values - coarsest first

   ,num_x     = imdi         &  ! Number of new coarse boxes x direction
                                ! for first new resolution
   ,num_y     = imdi         &  ! Number of new coarse boxes y direction
                                ! for first new resolution
   ,ntimes     = imdi        &  ! Number of times to be read in from input
                                ! files
   ,in_cols    = imdi        &  ! Number of columns in original fine grid
   ,in_rows    = imdi        &  ! Number of rows in original fine grid
   ,model_levels = imdi      &  ! Number of model levels in simulation
   ,bl_levels    = imdi      &  ! Number of BL levels in simulation
   ,mlevs      = imdi        &  ! Number of model levels in output being read
   ,num_ff     = imdi        &  ! Number of fieldsfiles to be read
   ,iprint      = imdi          ! 0 - reduced printed output
                                ! 1 - extra printed outpt to help check program

LOGICAL ::                   &
    l_bcu     = .FALSE.      &  ! .true. for BCU - buoyant cloudy updraughts
                                !                  (w'>0)
   ,l_wg1     = .FALSE.      &  ! .true. for wg1 - strong buoyant cloudy
                                !                  updraughts
   ,l_acc     = .FALSE.      &  ! .true. for ACC - all cloudy points
   ,l_acu     = .FALSE.      &  ! .true. for ACU - all cloudy updraughts (w'>0)
   ,l_ppd     = .FALSE.      &  ! .true. for PPD - precipitating downdraughts
                                !                  (w' <0)
   ,l_nbd     = .FALSE.      &  ! .true. for NBD - negatively buoyant
                                !            precipitating downdraughts(w' <0)
   ,l_nid     = .FALSE.      &  ! .true. for NID - negatively buoyant (+ice)
                                !                  precipitating downdraughts
   ,l_adu     = .FALSE.      &  ! .true. for ADU - Dry air moving upwards
   ,l_acw     = .FALSE.      &  ! .true. for ACW - Cloudy updraught (w>0)
   ,l_bcw     = .FALSE.      &  ! .true. for BCW - buoyant cloudy updraught
                                !                  (w>0)
   ,l_ucu     = .FALSE.      &  ! .true. for UCU - unstable cloudy updraughts
                                !               (w>0 & -dln(thetaesat)/dz >0)
   ,l_ppw     = .FALSE.      &  ! .true. for PPW precipitating downdraughts
                                !               (w<0)

   ,l_nbw     = .FALSE.      &  ! .true. for NBD - negatively buoyant
                                !            precipitating downdraughts (w<0)
   ,l_bcu_mask = .FALSE.     &  ! .true. for w at points which are buoyant
                                !     on (from analysis on first res)
   ,l_plume_size = .FALSE.   &  ! .true. for plume size analysis on coarse grid
   ,l_qgraup  = .FALSE.      &  ! .true. for  graupel in run
   ,l_pp = .FALSE.           &  ! .true. for pp rather than fieldsfile input
   ,l_whole_grid = .FALSE.   &  ! .true. for whole grid to be used suitable for
                                !        bicyclic crm style runs
   ,l_no_orog = .FALSE.      &  ! .true.  if no orography ancillary
   ,l_all_sea = .FALSE.      &  ! .true.  sea grid no land sea mask ancillary
   ,l_ENDGame = .FALSE.      &  ! .true. for ENDGame input grids
   ,l_class_col = .FALSE.    &  ! .true. work out column classification
                                !        based on surface precipitation rates
   ,l_sect30 = .FALSE.       &  ! .true. if section 30 diagnostics in input  
   ,l_cape   = .FALSE.       &  ! .true. if require CAPE & CIN calculations
   ,l_pcape  = .FALSE.          ! .true. if require PCAPE etc
                                
REAL  ::                     &
  adx         = rmdi         &  ! grid dx for idealised grid where no mask
 ,ady         = rmdi         &  ! grid dy
 ,azx         = rmdi         &  ! grid zx
 ,azy         = rmdi         &  ! grid zy
 ,apseudo_lat = rmdi         &  ! psuedo lat
 ,apseudo_lon = rmdi            ! psuedo lon

INTEGER ::                   &
  start_date(6)    = imdi    &  ! start date for input data - year
 ,modtstep_date(6) = imdi    &  ! model time step for input data - year
 ,datastep_date(6) = imdi       ! data time step for input data - year

!------------------------------------------------------------------------------
! Calculated from namelist input
!------------------------------------------------------------------------------
REAL  ::                     &
  timestep                      ! model time step (s) needed for dPCAPE/dt
                                ! calculations

!------------------------------------------------------------------------------
! Define namelist &Run_coarse_grid read in from control file.
!------------------------------------------------------------------------------

NAMELIST/run_coarse_grid/                                              &

! Integer
 nx_start, ny_start, nres, new_res, num_x, num_y, model_levels,        &
 bl_levels, mlevs, ntimes, num_ff, in_cols, in_rows, iprint,           &
! Logical
 l_bcu, l_wg1, l_acc, l_acu, l_ppd, l_nbd, l_nid, l_adu, l_acw, l_bcw, &
 l_ucu, l_ppw, l_nbw, l_bcu_mask, l_plume_size,                        &
 l_pp, l_qgraup, l_whole_grid, l_no_orog, l_all_sea, l_ENDGame,        &
 l_class_col, l_sect30, l_cape, l_pcape,                               &

! Real
  adx, ady, azy, azx, apseudo_lat, apseudo_lon

!------------------------------------------------------------------------------
! Define namelist &Run_date_info read in from control file.
!------------------------------------------------------------------------------

NAMELIST/run_date_info/                                                      &

! Integer
 start_date, modtstep_date, datastep_date

!------------------------------------------------------------------------------
! Control variables not set by namelist but used in code
!------------------------------------------------------------------------------

INTEGER, PARAMETER :: height_gen_method=2   ! choice of model heights

!------------------------------------------------------------------------------
! Fields required from input fields files
!------------------------------------------------------------------------------
INTEGER, PARAMETER :: num_want = 38   ! Number of fields required with sect30

INTEGER, PARAMETER :: num_want_no30 = 31  ! Number of fields required without 
                                          ! sect30

INTEGER, PARAMETER ::                                        &
  stash_list(38) = (/     2,     3,     4,    10,    12,     &
                         24,    25,   150,   254,   272,     &
                        273,   408,   409,  1181,  2181,     &
                       3184,  3217,  3234,  4181,  4182,     &
                       4183,  4184,  4203,  4204,  9181,     &
                       9182,  9183, 12181, 12182, 12183,     &
                      12184, 30181, 30182, 30183, 30184,     &
                      30188, 30189, 30190 /)                 &
 ,date_typ(38)   = (/     1,     1,     1,     1,     1,     &
                          1,     1,     1,     1,     1,     &
                          1,     1,     1,     2,     2,     &
                          2,     2,     2,     2,     2,     &
                          2,     2,     2,     2,     2,     &
                          2,     2,     2,     2,     2,     &
                          2,     2,     2,     2,     2,     &
                          2,     2,     2 /)                 &
 ,proc_list(38)   = (/    1,     2,     0,     0,     0,     &
                          0,     0,     0,     0,     0,     &
                          0,     3,     0,     0,     0,     &
                          0,     0,     0,     0,     0,     &
                          0,     0,     0,     0,     0,     &
                          0,     0,     0,     0,     0,     &
                          0,     0,     0,     0,     0,     &
                          0,     0,     0  /)

INTEGER ::                                                   &
  lev_list(38)   = (/    80,    80,    80,    80,    80,     &
                          1,     1,    80,    80,    80,     &
                         80,    80,     1,    80,    80,     &
                         80,     1,     1,    80,    80,     &
                         80,    80,     1,     1,    80,     &
                         80,    80,    80,    80,    80,     &
                         80,    80,    80,    80,    80,     &
                         80,    80,    80 /)
!------------------------------------------------------------------------------
! Info controlling pdfs of plumes - currently set here and not through
! a namelist
!------------------------------------------------------------------------------
INTEGER, PARAMETER :: nbins_diam = 30   ! Number of bins for diameters
INTEGER, PARAMETER :: nbins_size = 21   ! Number of bins for sizes
INTEGER, PARAMETER :: nbins_fract = 11  ! Number of bins for fractions
INTEGER, PARAMETER :: nbins_dxdy = 10   ! Number of bins for dxdy

!  diameters (30)
!  Every 1 up to 20, then every 5. gridpoints until 60.

REAL, PARAMETER ::  diam_bin(nbins_diam+1)                                   &
                                =(/0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,  &
                                  10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,   &
                                  20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,70.0,100.0/)

! Sizes in number of grid points
! 1., 1<=4, 4<=10, 20<=30, 30<=40, 40<=50, 50<=60, 60<=70, 70<=80, 80<=90,
! 90<=100, 100<=200, 300<=400, 300<=400, 400<=500, 500<=600, 600<=700,
! 700<=800, 800<=900, 900<=1000, >1000

REAL, PARAMETER :: &
  size_bin(nbins_size+1)=(/0.0,1.0,4.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,  &
                         80.0,90.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,   &
                         900.0,1000.0  /)


! Fraction
! <0.1, 0.1=<0.2, 0.2=<0.3, 0.3=<0.4, 0.4=<0.5,
! 0.5=<0.6, 0.6=<0.7, 0.7=<0.8, 0.8=<0.9, 0.9=<1.0,
! =1.0
REAL, PARAMETER ::         &
fract_bin(nbins_fract)=(/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)

!------------------------------------------------------------------------------
CONTAINS

SUBROUTINE print_nlist_run_coarse_grid()

USE umPrintMgr

IMPLICIT NONE

CALL umPrint('Contents of namelist run_coarse_grid',src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I10)') ' nx_start = ',nx_start
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I10)') ' ny_start = ',ny_start
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I10)') ' nres = ',nres
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I10)') ' new_res = ',new_res(10)
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I10)') ' num_x = ',num_x
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I10)') ' num_y = ',num_y
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I10)') ' ntimes = ',ntimes
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I10)') ' in_cols = ',in_cols
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I10)') ' in_rows = ',in_rows
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I10)') ' model_levels = ',model_levels
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,I10)') ' bl_levels = ',bl_levels
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,I10)') ' mlevs = ',mlevs
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,I10)') ' num_ff = ',num_ff
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,I10)') ' iprint = ',iprint
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,L2)') ' l_bcu = ',l_bcu
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_wg1 = ',l_wg1
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_acc = ',l_acc
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_acu = ',l_acu
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_ppd = ',l_ppd
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_nbd = ',l_nbd
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_nid = ',l_nid
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_adu = ',l_adu
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_acw = ',l_acw
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_bcw = ',l_bcw
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_ucu = ',l_ucu
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_ppw = ',l_ppw
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_nbw = ',l_nbw
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_bcu_mask = ',l_bcu_mask
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_plume_size = ',l_plume_size
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_qgraup = ',l_qgraup
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_pp = ',l_pp
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_whole_grid = ',l_whole_grid
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_no_orog = ',l_no_orog
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_all_sea = ',l_all_sea
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_ENDGame = ',l_ENDGame
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_class_col = ',l_class_col
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_cape = ',l_cape
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,L2)') ' l_pcape = ',l_pcape
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,F28.8)') ' adx = ',adx
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,F28.8)') ' ady = ',ady
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,F28.8)') ' azx = ',azx
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,F28.8)') ' azy = ',azy
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,F28.8)') ' apseudo_lat = ',apseudo_lat
CALL umPrint(umMessage,src='crmstyle_cntl_mod')
WRITE(umMessage,'(A,F28.8)') ' apseudo_lon = ',apseudo_lon
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

END SUBROUTINE print_nlist_run_coarse_grid

!------------------------------------------------------------------------------
SUBROUTINE print_nlist_run_date_info()


USE umPrintMgr

IMPLICIT NONE

CALL umPrint('Contents of namelist run_date_info',src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I5,5I3)') ' start_date    = ',start_date
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I5,5I3)') ' datastep_date = ',datastep_date
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

WRITE(umMessage,'(A,I5,5I3)') ' modtstep_date = ',modtstep_date
CALL umPrint(umMessage,src='crmstyle_cntl_mod')

END SUBROUTINE print_nlist_run_date_info
!------------------------------------------------------------------------------
END MODULE crmstyle_cntl_mod
