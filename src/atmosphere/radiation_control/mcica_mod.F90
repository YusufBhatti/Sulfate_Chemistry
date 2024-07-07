! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Global data module for variables used in McICA scheme
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------

MODULE mcica_mod

USE ereport_mod, ONLY: ereport
USE um_parcore, ONLY: mype
USE umprintmgr, ONLY: newline
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
SAVE


!     Variables Required for using MCICA scheme.
INTEGER :: ipph
!     Plane-parallel homogeneous flag

INTEGER :: ioverlap
!     Overlap flag

INTEGER :: tot_subcol_gen=64
!     number of sub-columns to generate

INTEGER,  ALLOCATABLE :: sw_subcol_k(:,:)
!     Number of subcolumns for each k_term in each band.

INTEGER,  ALLOCATABLE :: lw_subcol_k(:,:)
!     Number of subcolumns for each k_term in each band.

REAL, ALLOCATABLE  ::  frac_cloudy_full(:)
!     fraction of the profile which is cloudy, full array

INTEGER, ALLOCATABLE  ::  ncldy(:)
!     number of cloudy sub-columns in each profile

REAL, ALLOCATABLE  ::  clw_sub_full(:,:,:)
!     Value of cloud liquid water content for each sub-column

REAL, ALLOCATABLE  ::  cic_sub_full(:,:,:)
!     Value of cloud ice water content for each sub-column

INTEGER ::  subcol_need=59
!     Number of cloudy sub-columns required (i.e. MAX of SW and LW)
INTEGER ::  subcol_need_single
!     Number of cloudy sub-columns required for single sampling
INTEGER ::  subcol_need_optimal
!     Number of cloudy sub-columns required for optimal sampling

! Order of sub-columns points to either SW or LW. (Sub-columns are
! rearranged so that each sub-column is equivalently as important
! in the LW as in the SW.)
INTEGER, ALLOCATABLE ::  sw_subcol_reorder(:)
!     SW order of sub-columns

INTEGER, ALLOCATABLE ::  lw_subcol_reorder(:)
!     LW order of sub-columns

INTEGER, ALLOCATABLE ::  lw_subcol_reorder_single(:)
!     LW order of sub-columns (for single sampling)

INTEGER, ALLOCATABLE ::  lw_subcol_reorder_optimal(:)
!     LW order of sub-columns (for optimal sampling)

INTEGER, PARAMETER :: ip_mcica_full_sampling = 0
!     Each k-term "sees" every sub-column

INTEGER, PARAMETER :: ip_mcica_single_sampling = 1
!     Each k-term "sees" a single sub-column

INTEGER, PARAMETER :: ip_mcica_optimal_sampling = 2
!     Each k-term "sees" an optimal number of sub-columns

REAL, ALLOCATABLE ::  rand_seed_x(:, :)
!     global array of random numbers for first level in cloud generator

REAL, ALLOCATABLE ::  rand_seed_y(:, :)
!     global array of random numbers for first level in cloud generator

REAL, PARAMETER :: cut  = 0.001
!     Cutoff for minimum cloud amount in a layer

INTEGER :: n1, n2
!     Dimensions of xcw array:
!       Cumulative probability (n1)
!       Relative standard deviation (n2)

REAL, ALLOCATABLE :: xcw(:,:)
!     Distribution of normalised condensate amount as a function of
!     cumulative probability and relative standard deviation.

REAL, ALLOCATABLE :: cloud_inhom_param_full(:,:)
!     Scaling factor applied to water content to represent inhomogeneity
!     calculated from FSD parametrisation


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MCICA_MOD'

CONTAINS

SUBROUTINE read_mcica_data(mcica_data)

USE spec_sw_lw, ONLY: sw_spectrum, lw_spectrum
USE ereport_mod, ONLY: ereport
USE mpl,ONLY: mpl_packed, mpl_integer, mpl_real
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Intent IN arguments

CHARACTER (LEN=200), INTENT(IN) :: mcica_data
!     Path to McICA data file


! Local variables

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'READ_MCICA_DATA'
INTEGER, PARAMETER :: iu_mcd = 80
INTEGER            :: icode
CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=80) :: line
CHARACTER (LEN=errormessagelength) :: iomessage
INTEGER            :: band, k, subcol

INTEGER              :: my_comm
INTEGER              :: icode2
INTEGER              :: dimensions(10)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER              :: posn    
INTEGER              :: count_int
INTEGER              :: count_real
INTEGER              :: k2
INTEGER              :: k3
CHARACTER, ALLOCATABLE :: buffer(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! read in all information from file on task 0 and allocate appropriate arrays
IF (mype == 0) THEN 
  OPEN(UNIT=iu_mcd, FILE=mcica_data, IOSTAT=icode, STATUS='OLD',      &
       ACTION='READ',IOMSG=iomessage)
  IF (icode /= 0) THEN
    cmessage = 'McICA data file could not be opened.'//          newline//&
               'IoMsg: '//TRIM(iomessage)
    CALL ereport(routinename, icode, cmessage)
  END IF

  DO
    !     Read header until data block is found
    READ(iu_mcd, '(a80)', IOSTAT=icode, IOMSG=iomessage) line
    IF (line(1:5) == '*DATA') EXIT
    IF (icode /= 0) THEN
      cmessage = 'No *DATA block present in McICA data file: '// newline//&
                 'IoMsg: '//TRIM(iomessage)
      CALL ereport(routinename, icode, cmessage)
    END IF
  END DO

  DO
    !     Read variables from data block
    READ(iu_mcd, '(a80)', IOSTAT=icode) line
    IF (line(1:4) == '*END') EXIT

    SELECT CASE (line)

    CASE ('tot_subcol_gen')
      READ(iu_mcd, *, IOSTAT=icode) tot_subcol_gen
    CASE ('subcol_need_single')
      READ(iu_mcd, *, IOSTAT=icode) subcol_need_single
    CASE ('subcol_need_optimal')
      READ(iu_mcd, *, IOSTAT=icode) subcol_need_optimal
    CASE ('ipph')
      READ(iu_mcd, *, IOSTAT=icode) ipph

    CASE ('lw_subcol_reorder_single')
      ALLOCATE(lw_subcol_reorder_single(subcol_need_single),          &
        stat=icode)
      IF (icode /= 0) THEN
        cmessage = 'Cannot allocate array: lw_subcol_reorder_single'
        CALL ereport(routinename, icode, cmessage)
      END IF
      READ(iu_mcd, *, IOSTAT=icode) lw_subcol_reorder_single

    CASE ('lw_subcol_reorder_optimal')
      ALLOCATE(lw_subcol_reorder_optimal(subcol_need_optimal),        &
        stat=icode)
      IF (icode /= 0) THEN
        cmessage = 'Cannot allocate array: lw_subcol_reorder_optimal'
        CALL ereport(routinename, icode, cmessage)
      END IF
      READ(iu_mcd, *, IOSTAT=icode) lw_subcol_reorder_optimal

    CASE ('sw_subcol_k')
      ALLOCATE(sw_subcol_k(sw_spectrum(1)%dim%nd_band,                &
                         sw_spectrum(1)%dim%nd_k_term), stat=icode)
      IF (icode /= 0) THEN
        cmessage = 'Cannot allocate array: sw_subcol_k'
        CALL ereport(routinename, icode, cmessage)
      END IF
      sw_subcol_k=1
      DO
        READ(iu_mcd, '(3i4)', IOSTAT=icode, IOMSG=iomessage) band, k, subcol
        IF (band == -99) EXIT
        sw_subcol_k(band,k)=subcol
        IF (icode /= 0) THEN
          cmessage = 'Error reading data for sw_subcol_k:'// newline//&
                     'IoMsg: '//TRIM(iomessage)
          CALL ereport(routinename, icode, cmessage)
        END IF
      END DO
 
    CASE ('lw_subcol_k')
      ALLOCATE(lw_subcol_k(lw_spectrum(1)%dim%nd_band,                &
                         lw_spectrum(1)%dim%nd_k_term), stat=icode)
      IF (icode /= 0) THEN
        cmessage = 'Cannot allocate array: lw_subcol_k'
        CALL ereport(routinename, icode, cmessage)
      END IF
      lw_subcol_k=1
      DO
        READ(iu_mcd, '(3i4)', IOSTAT=icode, IOMSG=iomessage) band,k, subcol
        IF (band == -99) EXIT
        lw_subcol_k(band,k)=subcol
        IF (icode /= 0) THEN
          cmessage = 'Error reading data for lw_subcol_k: '// newline//&
                     'IoMsg: '//TRIM(iomessage)
          CALL ereport(routinename, icode, cmessage)
        END IF
      END DO

    CASE ('n1')
      READ(iu_mcd, *, IOSTAT=icode) n1
    CASE ('n2')
      READ(iu_mcd, *, IOSTAT=icode) n2
    CASE ('xcw')
      ALLOCATE(xcw(n1,n2), stat=icode)
      IF (icode /= 0) THEN
        cmessage = 'Cannot allocate array: xcw'
        CALL ereport(routinename, icode, cmessage)
      END IF
      READ(iu_mcd, *, IOSTAT=icode) xcw
  
    END SELECT
 
    IF (icode /= 0) THEN
      cmessage = 'Error reading data from McICA data file'
      CALL ereport(routinename, icode, cmessage)
    END IF
  END DO

  CLOSE(iu_mcd)

END IF  !mype==0

CALL gc_get_communicator(my_comm, icode2)

! broadcast all integers and array size information from task 0 to all others
! to allow array allocation

IF (mype == 0) THEN
  dimensions(1)  = SIZE(sw_subcol_k, 1)
  dimensions(2)  = SIZE(sw_subcol_k, 2)
  dimensions(3)  = SIZE(lw_subcol_k, 1)
  dimensions(4)  = SIZE(lw_subcol_k, 2)
  dimensions(5)  = subcol_need_single
  dimensions(6)  = subcol_need_optimal
  dimensions(7)  = n1
  dimensions(8)  = n2
  dimensions(9)  = ipph
  dimensions(10) = tot_subcol_gen
END IF

CALL mpl_bcast(dimensions, 10, mpl_integer, 0, my_comm, icode2)

IF (mype /= 0) THEN
  ALLOCATE ( sw_subcol_k(dimensions(1), dimensions(2)) )
  ALLOCATE ( lw_subcol_k(dimensions(3), dimensions(4)) )
  ALLOCATE ( lw_subcol_reorder_single(dimensions(5)) )
  ALLOCATE ( lw_subcol_reorder_optimal(dimensions(6)) )
  ALLOCATE ( xcw(dimensions(7), dimensions(8)) )
  subcol_need_single  = dimensions(5)
  subcol_need_optimal = dimensions(6)
  n1             = dimensions(7)
  n2             = dimensions(8)
  ipph           = dimensions(9)
  tot_subcol_gen = dimensions(10)
END IF

! now size of arrays known can allocate buffer size required to transmit data
! note extra size for safety

count_int = dimensions(1)*dimensions(2) + dimensions(3)*dimensions(4) &
           + dimensions(5) + dimensions(6)
count_real = dimensions(7)*dimensions(8)   

CALL mpl_pack_size(count_int, mpl_integer, my_comm, k2, icode2)
CALL mpl_pack_size(count_real, mpl_real, my_comm, k3, icode2)
k=k2+k3
k= INT( (11*k) /10 )    !! add 10% to buffer for safety

ALLOCATE(buffer(k))
posn=0

IF ( mype == 0 ) THEN
  CALL mpl_pack(sw_subcol_k, dimensions(1)*dimensions(2),           &
               mpl_integer, buffer, k, posn, my_comm, icode2)
  CALL mpl_pack(lw_subcol_k, dimensions(3)*dimensions(4),           &
               mpl_integer, buffer, k, posn, my_comm, icode2)
  CALL mpl_pack(lw_subcol_reorder_single, dimensions(5),            &
               mpl_integer, buffer, k, posn, my_comm, icode2)
  CALL mpl_pack(lw_subcol_reorder_optimal, dimensions(6),           &
               mpl_integer, buffer, k, posn, my_comm, icode2)
  CALL mpl_pack(xcw, dimensions(7)*dimensions(8),                   &
               mpl_real, buffer, k, posn, my_comm, icode2)
END IF  ! mype==0

! broadcast buffer
CALL mpl_bcast(buffer, k, mpl_packed, 0, my_comm, icode2)

! unpack data on pe/=0
IF ( mype /= 0 ) THEN
 
  CALL mpl_unpack(buffer, k, posn, sw_subcol_k,                 &
     dimensions(1)*dimensions(2),  mpl_integer, my_comm, icode2)
  CALL mpl_unpack(buffer, k, posn, lw_subcol_k,                 &
     dimensions(3)*dimensions(4),  mpl_integer, my_comm, icode2)
  CALL mpl_unpack(buffer, k, posn, lw_subcol_reorder_single,    &
     dimensions(5),  mpl_integer, my_comm, icode2)
  CALL mpl_unpack(buffer, k, posn, lw_subcol_reorder_optimal,   &
     dimensions(6),  mpl_integer, my_comm, icode2)
  CALL mpl_unpack(buffer, k, posn, xcw,                         &
     dimensions(7)*dimensions(8),  mpl_real, my_comm, icode2)
END IF  ! mype==0

DEALLOCATE(buffer)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE read_mcica_data

END MODULE mcica_mod
