! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stash

MODULE cstash_mod

USE version_mod, ONLY: ndiagpm, nproftp, ntimep, nprofdp, NTimSerP,  &
                       nlevp, npslevp,npslistp,nprofup,nprofpp

USE profilename_length_mod, ONLY: profilename_length,packagename_length

IMPLICIT NONE

! Description:
!  Contains variables and arrays involved in STASH request specification
!  and STASH processing in the UM. Includes namelist STASH.

! Global scalars:
INTEGER ::   ndiag   ! No. of diagnostics
INTEGER ::   ntprof  ! No. of time profiles
INTEGER ::   nseries ! No. of stash time series
INTEGER ::   ndprof  ! No. of domain profiles
INTEGER ::   nuprof  ! No. of useage profiles
INTEGER ::   npprof  ! No. of excluded package profiles

! Global dynamic arrays:

!   STASH specification table (JSTASH file):
!   NDIAGPM set in VERSION_MOD module
INTEGER ::   modl_b(ndiagpm)  ! Internal model no.
INTEGER ::   isec_b(ndiagpm)  ! Section
INTEGER ::   item_b(ndiagpm)  ! Item
INTEGER ::   itim_b(ndiagpm)  ! Time profile number
INTEGER ::   idom_b(ndiagpm)  ! Domain profile number
INTEGER ::   iuse_b(ndiagpm)  ! Useage profile number

!   Time profile information:

CHARACTER(LEN=profilename_length) :: timpro(nproftp)   ! Name of time profile
INTEGER ::   ityp_t(nproftp)         ! Type of profile
INTEGER ::   intv_t(nproftp)         ! Time Interval
INTEGER ::   unt1_t(nproftp)         ! Units for time interval
INTEGER ::   isam_t(nproftp)         ! Sampling period
INTEGER ::   unt2_t(nproftp)         ! Units for sampling period
INTEGER ::   iopt_t(nproftp)         ! Output option
INTEGER ::   istr_t(nproftp)         ! Output Start time
INTEGER ::   iend_t(nproftp)         ! Output End time
INTEGER ::   isdt_t(6, nproftp)      ! Output Start date
INTEGER ::   iedt_t(6, nproftp)      ! Output End date
INTEGER ::   ifre_t(nproftp)         ! Output frequency
INTEGER ::   ioff_t(nproftp)         ! Offset for sampling
INTEGER ::   unt3_t(nproftp)         ! Units for output times
INTEGER ::   itim_t(nproftp)         ! No. of times in times table
INTEGER ::   iser_t(ntimep ,nproftp) ! Times table (with units)
INTEGER ::   modl_t(nproftp)         ! Indicates internal model
                                     !  for each times table
LOGICAL ::   lts0_t(nproftp)         ! Logical for PP output at timestep zero 

!   Domain profile information:
CHARACTER(LEN=profilename_length) :: dompro  (nprofdp) ! Name of domain profile
INTEGER ::  iopl_d  (nprofdp)           ! Levels option
INTEGER ::  levb_d  (nprofdp)           ! Bottom level
INTEGER ::  levt_d  (nprofdp)           ! Top level
INTEGER ::  iopa_d  (nprofdp)           ! Area option
INTEGER ::  inth_d  (nprofdp)           ! North boundary
INTEGER ::  isth_d  (nprofdp)           ! South boundary
INTEGER ::  iest_d  (nprofdp)           ! East boundary
INTEGER ::  iwst_d  (nprofdp)           ! West boundary
INTEGER ::  imsk_d  (nprofdp)           ! Mask type
INTEGER ::  imn_d   (nprofdp)           ! Meaning option
INTEGER ::  iwt_d   (nprofdp)           ! Weighting option
LOGICAL ::  ts_d    (nprofdp)           ! Time series profile
INTEGER ::  ig_ts
INTEGER ::  i1_ts
INTEGER ::  i51_ts
INTEGER ::  blim_ts (NTimSerP)
INTEGER ::  tlim_ts (NTimSerP)
REAL ::     blimr_ts(NTimSerP)
REAL ::     tlimr_ts(NTimSerP)
INTEGER ::  nlim_ts (NTimSerP)
INTEGER ::  slim_ts (NTimSerP)
INTEGER ::  elim_ts (NTimSerP)
INTEGER ::  wlim_ts (NTimSerP)
INTEGER ::  ilev_d  (nprofdp)           ! Output levels code
INTEGER ::  levlst_d(nlevp   ,nprofdp ) ! Levels list
REAL ::    rlevlst_d(nlevp   ,nprofdp ) ! Levels list
INTEGER ::  plt_d   (nprofdp)
INTEGER ::  pllen_d (nprofdp)
INTEGER ::  plpos_d (nprofdp)
INTEGER ::  pslist_d(npslevp ,npslistp)
INTEGER ::  npslists
EQUIVALENCE        (rlevlst_d,levlst_d)

! Useage information:

CHARACTER(LEN=profilename_length) :: usepro(nprofup) ! Name of useage profile
INTEGER ::  locn_u(nprofup)   ! Storage location of profile
INTEGER ::  iunt_u(nprofup)   ! Unit no.
LOGICAL ::  lnetcdf_u(nprofup)! Is output format netCDF?

!Package information
CHARACTER(LEN=packagename_length) :: pckpro(nprofpp)

! Information from ppxref file:

INTEGER ::   model_st       ! Internal model number
INTEGER ::   ispace         ! Space code
INTEGER ::   itima          ! Time availability code
INTEGER ::   igp            ! Grid of data code
INTEGER ::   ilev           ! Level type code
INTEGER ::   ibot           ! First level code
INTEGER ::   itop           ! Last level code
INTEGER ::   iflag          ! Level compression flag
INTEGER ::   iopn(6)        ! Sectional option code
INTEGER ::   vmsk           ! Integer equiv of bin vers mask
INTEGER ::   ipseudo        ! Pseudo dimension type
INTEGER ::   ipfirst        ! First pseudo dim code
INTEGER ::   iplast         ! Last pseudo dim code
INTEGER ::   ptr_prog       ! Section zero point back
INTEGER ::   halo_type      ! Type of halo the field has

! Available time units

INTEGER, PARAMETER ::  stsh_timesteps = 1
INTEGER, PARAMETER ::  stsh_hours     = 2
INTEGER, PARAMETER ::  stsh_days      = 3
INTEGER, PARAMETER ::  stsh_dumps     = 4
INTEGER, PARAMETER ::  stsh_minutes   = 5
INTEGER, PARAMETER ::  stsh_seconds   = 6

END MODULE cstash_mod
