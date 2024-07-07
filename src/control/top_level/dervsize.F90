! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Program: DERVSIZE -------------------------------------------------
!
!    Purpose: Calculate extra sizes required for dynamic allocation of
!             main memory in the model, derived from sizes passed by
!             READSIZE into the top level program UM_SHELL. These are
!             local sizes for each pe, except where explicitly stated,
!             having previously called decomposition routines.
!
!    Programming standard: UM Doc Paper 3
!
!    External documentation: On-line UM document C1 - Dynamic allocation
!                            of primary fields
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE dervsize(                                              &
             icode,cmessage)
!
! ----------------------------------------------------------------------

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: atm_fields_bounds_init, tdims 
USE umPrintMgr
USE UM_ParVars
USE UM_ParCore,          ONLY: nproc
USE Control_Max_Sizes
USE Decomp_DB
USE lbc_mod
USE dust_parameters_mod, ONLY:               l_dust,              &
     l_dust_div1,       l_dust_div2,         l_dust_div3,         &
     l_dust_div4,       l_dust_div5,         l_dust_div6,         &
     l_dust_div1_lbc,   l_dust_div2_lbc,     l_dust_div3_lbc,     &
     l_dust_div4_lbc,   l_dust_div5_lbc,     l_dust_div6_lbc


USE run_aerosol_mod, ONLY: l_ocff_new, l_ocff_agd, l_ocff_cld,  &
                           l_soot_new, l_soot_agd, l_soot_cld,  &
                        l_bmass_new, l_bmass_agd, l_bmass_cld,  &
                  l_so4_aitken, l_so4_accu, l_so4_diss, l_so2,  &
                  l_dms, l_nh3, l_nitr_acc, l_nitr_diss,        &
                  l_so2_lbc, l_dms_lbc, l_so4_aitken_lbc,       &
                  l_so4_accu_lbc, l_so4_diss_lbc, l_nh3_lbc,    &
              l_soot_new_lbc, l_soot_agd_lbc, l_soot_cld_lbc,   &
              l_bmass_new_lbc, l_bmass_agd_lbc, l_bmass_cld_lbc,&
              l_ocff_new_lbc, l_ocff_agd_lbc, l_ocff_cld_lbc,   &
              l_nitr_acc_lbc, l_nitr_diss_lbc

USE mphys_inputs_mod, ONLY:                                        &
     l_mcr_qcf2,        l_mcr_qrain,         l_mcr_qgraup,        &
     l_mcr_qcf2_lbc,    l_mcr_qrain_lbc,     l_mcr_qgraup_lbc

USE cloud_inputs_mod, ONLY: l_pc2_lbc
USE murk_inputs_mod,  ONLY: l_murk, l_murk_lbc
USE cv_run_mod,       ONLY: l_3d_cca, l_ccrad

USE lbc_read_data_mod, ONLY: l_old_lbc_file, L_int_uvw_lbc

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    bl_levels, global_row_length, global_rows, intf_len2_coldepc,      &
    intf_len2_levdepc, intf_len2_rowdepc, intf_lookupsa, model_levels, &
    n_cca_lev, n_rows, ozone_levels, row_length, rows,                 &
    theta_field_size, theta_halo_size, theta_off_size, tr_levels,      &
    u_field_size, u_halo_size, u_off_size, v_field_size, v_halo_size,  &
    v_off_size

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_bi_cyclic_lam

IMPLICIT NONE
!
!  Argument list
!
INTEGER :: icode             ! OUT - Return code
CHARACTER(LEN=errormessagelength) :: cmessage     ! OUT - Error message

! For STASH sizes
!!! #include "typstsz.h"
!  Local variables
INTEGER ::                                                        &
  iside                                                           &
               ! loop counter for LBC sides
, ifld                                                            &
               ! loop counter for field types
, ihalo                                                           &
               ! loop counter for halo types
, iproc                                                           &
               ! loop counter for processors
, irim                                                            &
               ! loop counter for rim types
, info                                                            &
               ! return code from GCOM
, lbc_row_len                                                     &
               ! length of row of LBC
, lbc_nrows                                                       &
               ! number of rows in LBC
, num_optional_lbcs_in  ! no. of optional lbc fields in input

! ----------------------------------------------------------------------

INTEGER :: nohalo_IMT,nohalo_JMT,glob_IMT,glob_JMT

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DERVSIZE'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
icode=0

!     Initialise number of optional input lbcs to zero
num_optional_lbcs_in  = 0
!
!    Atmosphere Boundary Datasets.
!    2nd dimension of Level Dependent Constants.
intf_len2_levdepc=4
!    2nd dimension of Row/Col Dependent Constants.
intf_len2_rowdepc=2
intf_len2_coldepc=2

! Sizes applicable to all resolutions
theta_field_size=row_length*rows
! One less V row on the Northern most processors
IF (at_extremity(PNorth) .AND.                                    &
         model_type /= mt_bi_cyclic_lam) THEN
  n_rows=rows+1
ELSE
  n_rows=rows
END IF
u_field_size=theta_field_size
v_field_size=row_length*n_rows
theta_off_size   = (row_length + 2*offx)   * (rows   + 2*offy)
theta_halo_size  = (row_length + 2*halo_i) * (rows   + 2*halo_j)
u_off_size       = (row_length + 2*offx)   * (rows   + 2*offy)
u_halo_size      = (row_length + 2*halo_i) * (rows   + 2*halo_j)
v_off_size       = (row_length + 2*offx)   * (n_rows + 2*offy)
v_halo_size      = (row_length + 2*halo_i) * (n_rows + 2*halo_j)

!     Grid bounds settings
CALL atm_fields_bounds_init(offx,offy,halo_i,halo_j,              &
            row_length,rows,n_rows,                               &
            tr_levels,bl_levels,ozone_levels)

!     No of levels for Convective Cloud Amount.
IF (l_3d_cca .OR. l_ccrad) THEN
  ! This needs to be the number of wet levels from level 1 - even for EG.
  n_cca_lev = tdims%k_end
ELSE
  n_cca_lev = 1
END IF
IF (PrintStatus >= PrStatus_Normal) THEN
  WRITE(umMessage,*)                                                      &
  'DERVSIZE: Number of levels for convective clouds is ',         &
  n_cca_lev
  CALL umPrint(umMessage,src='dervsize')
END IF

! ----------------------------------------------------------------
! Count number of optional lateral boundary fields expected in
! input dependent on whether the _lbc logicals are true or false
! ----------------------------------------------------------------

! Additional microphysics variables (ice crystals, rain, graupel)
IF (L_mcr_qcf2_lbc) THEN  ! qcf2 lbcs active
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
IF (L_mcr_qrain_lbc) THEN  ! qrain lbcs active
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
IF (L_mcr_qgraup_lbc) THEN  ! qgraup lbcs active
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

! Cloud fractions for PC2 (3 fields: bulk, liquid and frozen)
IF (L_pc2_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 3
END IF

! Murk aerosol
IF (L_murk_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

! Dust
IF (L_dust_div1_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
IF (L_dust_div2_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
IF (L_dust_div3_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
IF (L_dust_div4_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
IF (L_dust_div5_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
IF (L_dust_div6_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

! SO2
IF (L_so2_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

! DMS
IF (L_dms_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

! SO4_Aitken
IF (L_so4_aitken_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
! SO4_ACCU
IF (L_so4_accu_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
! SO4_DISS
IF (L_so4_diss_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

! NH3
IF (L_NH3_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

! SOOT_NEW
IF (L_soot_new_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
! SOOT_AGD
IF (L_soot_agd_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
! SOOT_CLD
IF (L_soot_cld_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

! BMASS_NEW
IF (L_bmass_new_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
! BMASS_AGD
IF (L_bmass_agd_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
! BMASS_CLD
IF (L_bmass_cld_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

! OCFF_NEW
IF (L_ocff_new_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
! OCFF_AGD
IF (L_ocff_agd_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
! OCFF_CLD
IF (L_ocff_cld_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

! NITR_ACC
IF (L_nitr_acc_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF
! NITR_DISS
IF (L_nitr_diss_lbc) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 1
END IF

!Treat INPUT advected winds as optional
IF (.NOT. L_int_uvw_lbc .OR. L_old_lbc_file) THEN
  num_optional_lbcs_in = num_optional_lbcs_in + 3
END IF

! ----------------------------------------------------------------

a_len1_levdepc=model_levels+1
! We use the global values here
a_len1_rowdepc= global_rows + 1

a_len1_coldepc= global_row_length

!  Number of atmosphere model interface lookups
intf_lookupsa = 10

!
! Calculate sizes of Lateral Boundary Conditions Data:
! LENRIMA(fld_type,halo_type,rim_type)
!    size of a single level of LBC data for a given field type and
!    halo type and rim_type
! LBC_SIZEA(side,fld_type,halo_type,rim_type)
!    size of a single edge of LBC data for a given edge (North,
!    East,South,West) and a given field type and halo type and rim type
! LBC_STARTA(side,fld_type,halo_type,rim_type)
!    start address in LBC array of a given edge for a given
!    field type and halo type and rim_type
!

IF (.NOT. ASSOCIATED(g_lbc_starta)) &
    ALLOCATE(g_lbc_starta          &
    (4,nfld_max,nhalo_max,nrima_max,0:nproc_max-1))

IF (.NOT. ASSOCIATED(g_lenrima)) &
    ALLOCATE(g_lenrima          &
    (nfld_max,nhalo_max,nrima_max,0:nproc_max-1))

DO irim=1,Nrima_max
  DO ifld=1,Nfld_max        ! loop over field types
    DO ihalo=1,NHalo_max    ! loop over halo types

      lenrima(ifld,ihalo,irim)=0
      global_LENRIMA(ifld,ihalo,irim)=0

      DO iproc=0,nproc-1
        g_lenrima(ifld,ihalo,irim,iproc)=0
      END DO

      DO iside=1,4          ! loop over North,East,South,West

        DO iproc=0,nproc-1
          g_LBC_STARTA(iside,ifld,ihalo,irim,iproc)=0
        END DO
        ! First calculate the global_LENRIMA values - these are the sizes
        ! of the complete LBC as stored on disk before it is decomposed
        ! over processors.

        IF ((iside  ==  PNorth) .OR. (iside  ==  PSouth)) THEN
          ! North or South boundaries

          lbc_row_len=glsize(1,ifld)

          lbc_row_len=lbc_row_len+2*halosize(1,ihalo)

          lbc_nrows=halosize(2,ihalo)+rimwidtha(irim)


        ELSE
          ! East or West boundaries
          lbc_row_len=halosize(1,ihalo)+rimwidtha(irim)

          lbc_nrows=glsize(2,ifld)-2*rimwidtha(irim)

        END IF ! North/South or East/West boundary

        IF (rimwidtha(irim)  >   0) THEN

          global_LBC_SIZEA(iside,ifld,ihalo,irim)=                  &
            lbc_row_len*lbc_nrows

          IF (iside  ==  1) THEN
            global_LBC_STARTA(iside,ifld,ihalo,irim)=1
          ELSE
            global_LBC_STARTA(iside,ifld,ihalo,irim)=               &
              global_LBC_STARTA(iside-1,ifld,ihalo,irim)+           &
              global_LBC_SIZEA(iside-1,ifld,ihalo,irim)
          END IF

          global_LENRIMA(ifld,ihalo,irim)=                          &
            global_LENRIMA(ifld,ihalo,irim)+                        &
            global_LBC_SIZEA(iside,ifld,ihalo,irim)

        ELSE ! No LBCs if RIMWIDTH is  <=  0)

          global_LBC_SIZEA(iside,ifld,ihalo,irim)=0
          global_LBC_STARTA(iside,ifld,ihalo,irim)=1

        END IF ! IF (RIMWIDTHA  >   0)

        ! Now calculate local LENRIMA values, and the associated arrays
        ! LBC_SIZEA and LBC_STARTA

        IF (at_extremity(iside) .AND.                               &
            (rimwidtha(irim)  >   0)) THEN
          ! This processor is at the edge of the grid

          IF ((iside  ==  PNorth) .OR. (iside  ==  PSouth)) THEN
            ! North or South boundaries
            ! North/South boundaries can contain the corners of the LBCs

            ! Length of rows includes the halos.
            lbc_row_len=lasize(1,ifld,ihalo)

            ! And the number of rows is the size of the halo plus the rimwidth
            lbc_nrows=halosize(2,ihalo)+rimwidtha(irim)

          ELSE
            ! East or West boundaries

            ! Length of row is the size of the halo plus the rimwidth
            lbc_row_len=halosize(1,ihalo)+rimwidtha(irim)

            ! Number of rows is the number of rows on this PE minus any
            ! rows which are looked after by the North/South boundaries
            ! (ie. the corners).
            lbc_nrows=lasize(2,ifld,ihalo)
            IF (at_extremity(PNorth))                               &
              lbc_nrows=lbc_nrows-halosize(2,ihalo)-rimwidtha(irim)
            IF (at_extremity(PSouth))                               &
              lbc_nrows=lbc_nrows-halosize(2,ihalo)-rimwidtha(irim)

          END IF ! North/South or East/West boundary

          lbc_sizea(iside,ifld,ihalo,irim)=lbc_row_len*lbc_nrows

        ELSE
          ! This processor is not at the edge of the grid, or RIMWIDTH is
          ! zero (indicating no LBCs)
          lbc_sizea(iside,ifld,ihalo,irim)=0

        END IF

        ! LBC_STARTA contains the offset in the LBC array for each side
        ! (North,East,South,West) piece of data

        IF (iside  ==  1) THEN
          lbc_starta(iside,ifld,ihalo,irim)=1
        ELSE
          lbc_starta(iside,ifld,ihalo,irim)=                        &
            lbc_starta(iside-1,ifld,ihalo,irim)+                    &
            lbc_sizea(iside-1,ifld,ihalo,irim)
        END IF
        g_LBC_STARTA(iside,ifld,ihalo,irim,mype)=                   &
          lbc_starta(iside,ifld,ihalo,irim)

        lenrima(ifld,ihalo,irim)=lenrima(ifld,ihalo,irim)+           &
                                 lbc_sizea(iside,ifld,ihalo,irim)

        g_lenrima(ifld,ihalo,irim,mype)=                   &
          lenrima(ifld,ihalo,irim)

      END DO ! iside
    END DO ! ihalo
  END DO ! ifld
END DO ! irim

! Now do some comms so that g_LBC_STARTA is filled with the
! value of LBC_STARTA on every processor

CALL gc_imax(4*Nfld_max*NHalo_max*Nrima_max*nproc,nproc,info,     &
             g_LBC_STARTA)

! And the same for g_lenrima

CALL gc_imax(Nfld_max*NHalo_max*Nrima_max*nproc,nproc,info,     &
             g_lenrima)

! And set up a few other variables

      ! Includes one-off orography field at start
      ! Will be updated later, when we know the number of tracer LBCs
rim_lookupsa = 10 + num_optional_lbcs_in

!   LAM DERIVED SIZES END
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE dervsize
