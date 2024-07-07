! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   To put fields required by RADAER into stash.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE glomap_clim_fields_mod

USE atm_fields_real_mod,             ONLY: &
    cf_bulk,                               &
    cf_liquid,                             &
    exner_theta_levels,                    &
    gc_nd_nuc_sol,                         &
    gc_nuc_sol_su,                         &
    gc_nuc_sol_oc,                         &
    gc_nd_ait_sol,                         &
    gc_ait_sol_su,                         &
    gc_ait_sol_bc,                         &
    gc_ait_sol_oc,                         &
    gc_nd_acc_sol,                         &
    gc_acc_sol_su,                         &
    gc_acc_sol_bc,                         &
    gc_acc_sol_oc,                         &
    gc_acc_sol_ss,                         &
    gc_nd_cor_sol,                         &
    gc_cor_sol_su,                         &
    gc_cor_sol_bc,                         &
    gc_cor_sol_oc,                         &
    gc_cor_sol_ss,                         &
    gc_nd_ait_ins,                         &
    gc_ait_ins_bc,                         &
    gc_ait_ins_oc,                         &
    p_theta_levels,                        &
    q,                                     &
    qcf,                                   &
    theta

USE ereport_mod,                     ONLY: &
    ereport

USE errormessagelength_mod,          ONLY: &
    errormessagelength

USE ukca_mode_setup,                 ONLY: &
    component,                             &
    cp_su,                                 &
    cp_bc,                                 &
    cp_oc,                                 &
    cp_cl,                                 &
    cp_so,                                 &
    cp_du,                                 &
    mode_nuc_sol,                          &
    mode_ait_sol,                          &
    mode_acc_sol,                          &
    mode_cor_sol,                          &
    mode_ait_insol,                        &
    mode_acc_insol,                        &
    mode_cor_insol,                        &
    mfrac_0,                               &
    mm,                                    &
    mmid,                                  &
    mlo,                                   &
    mode,                                  &
    ncp,                                   &
    nmodes,                                &
    num_eps

USE umPrintMgr,                      ONLY: &
    umPrint,                               &
    umMessage,                             &
    PrintStatus,                           &
    PrStatus_Diag


IMPLICIT NONE
PRIVATE

! -----------------------------------------------------------------------------*
!                                                                              *
! Procedure:                                                                   *
!   1) Identify GLOMAP-mode climatology aerosol fields. Also identify          *
!       fields to be passed to RADAER (which need to be copied to stash)       *
!                                                                              *
!   2) Read in aerosol and other fields from D1                                *
!                                                                              *
!   3) CALL ukca_calc_drydiam and ukca_volume_mode routines to                 *
!       calculate the fields required by RADAER                                *
!                                                                              *
! Included subroutines:                                                        *
!         prepare_fields_for_radaer     (public)                               *
!         get_d1_fields                 (public)                               *
!         set_field                     (private)                              *
! -----------------------------------------------------------------------------*

REAL, ALLOCATABLE, PUBLIC :: mmr3d(:,:,:,:,:)
! Avg cpt mass mixing ratio of aerosol particle in mode (particle^-1)

REAL, ALLOCATABLE, PUBLIC :: nmr3d(:,:,:,:)
! Aerosol ptcl (number density/ air density) for mode (cm^-3)

REAL, ALLOCATABLE, PUBLIC :: md(:,:,:)
! Avg cpt mass of aerosol particle in mode (particle^-1) (n_points,nmodes,ncp)

REAL, ALLOCATABLE, PUBLIC :: nd(:,:)
! Aerosol ptcl number density for mode (cm^-3) (n_points,nmodes)

REAL, ALLOCATABLE, PUBLIC :: mdt(:,:)
! Avg tot mass of aerosol ptcl in mode (particle^-1)  (n_points,nmodes)

REAL, ALLOCATABLE, PUBLIC :: drydp(:,:)
! Median particle dry diameter for each mode (m)      (n_points,nmodes)

REAL, ALLOCATABLE :: wetdp(:,:)
! Avg wet diameter of size mode (m)                   [n_points,nmodes]

REAL, ALLOCATABLE :: rhopar(:,:)
! Particle density (kg/m^3) [includes H2O & insoluble cpts]

REAL, ALLOCATABLE, PUBLIC :: dvol(:,:)
! Median particle dry volume for each mode (m^3)      (n_points,nmodes)

REAL, ALLOCATABLE, PUBLIC :: aird(:)
! Air density number concentration (/cm3) (n_points)

REAL, ALLOCATABLE :: wvol(:,:)
! Avg wet volume of size mode (m3)  (n_points,nmodes)

REAL, ALLOCATABLE :: mdwat(:,:)
! Molecular concentration of water (molecules/particle) (n_points,nmodes)

REAL, ALLOCATABLE :: pvol_wat(:,:)
! Partial volume of water in each mode (m3)  (n_points,nmodes)

REAL, ALLOCATABLE :: pvol(:,:,:)
! Partial volumes of each component in each mode (m3)   (n_points,nmodes,ncp)

REAL, ALLOCATABLE :: qsatmr(:,:)
! Saturated specific humidity

REAL, ALLOCATABLE :: qsatmr_wat(:,:)
! Sat. mixing ratio with respect to liquid water irrespective of temperature

REAL, ALLOCATABLE :: t_theta_levels(:,:,:)
! Temperature on theta levels without halos

REAL, ALLOCATABLE, PRIVATE :: pr_theta_levels(:,:,:)
! Pressure at mid levels (3-D) (Pa)

REAL, ALLOCATABLE :: rel_humid_frac_clr(:,:,:)
! Clear sky relative humidity as a fraction

REAL, ALLOCATABLE, PUBLIC :: pmid_1d(:)
! Pressure at mid levels (Pa)

REAL, ALLOCATABLE :: q_1d(:)
! Specific humidity 1-D

REAL, ALLOCATABLE :: qcf_1d(:)
! qcf (kg/kg) 1-D

REAL, ALLOCATABLE :: q_clr_1d(:)
! Clear sky specific humidity 1-D

REAL, ALLOCATABLE :: qsatmr_1d(:)
! Saturated specific humidity 1-D

REAL, ALLOCATABLE :: qsatmr_wat_1d(:)
! Saturated specific humidity w.r.t. water 1-D

REAL, ALLOCATABLE :: rh_clr_1d(:)
! Clear sky relative humidity as a fraction 1-D

REAL, ALLOCATABLE :: rhcrit_1d(:)
! RH crit 1D

REAL, ALLOCATABLE :: cloud_liq_frac_1d(:)
! Cloud liquid fraction 1-D

REAL, ALLOCATABLE :: cloud_blk_frac_1d(:)
! Cloud bulk fraction 1-D

REAL, ALLOCATABLE, PUBLIC :: temp_1d(:)
! Temperature 1-D

CHARACTER(LEN=*), PARAMETER :: ModuleName='GLOMAP_CLIM_FIELDS_MOD'

PUBLIC :: prepare_fields_for_radaer, get_d1_fields

CONTAINS

! ##############################################################################

SUBROUTINE prepare_fields_for_radaer( ukca_radaer, stashwork54 )

! Read in mass mixing ratio and number from D1 to create the md and nd arrays,
!  THEN CALL GLOMAP routines to establish the fields required by RADAER.
!  Fill the section 54 arrays required by RADAER.

USE atm_fields_bounds_mod,         ONLY: &
    tdims

USE cloud_inputs_mod,              ONLY: &
    rhcrit

USE glomap_clim_identify_fields_mod, ONLY: &
    aerofields,                            &
    glomap_clim_identify_fields,           &
    n_stored_items

USE glomap_clim_option_mod,        ONLY: &
    l_glomap_clim_arg_act

USE lsp_subgrid_mod,               ONLY: &
    lsp_qclear

USE missing_data_mod,              ONLY: &
    imdi

USE nlsizes_namelist_mod,          ONLY: &
    row_length,                          &
    rows,                                &
    model_levels

USE parkind1,                      ONLY: &
    jprb,                                &
    jpim

USE qsat_mod,                      ONLY: &
    qsat,                                &
    qsat_wat_mix

USE stash_array_mod,               ONLY: &
    len_stlist,                          &
    num_stash_levels,                    &
    stlist,                              &
    stash_maxlen,                        &
    sf,                                  &
    si,                                  &
    stindex,                             &
    stash_levels

USE submodel_mod,                  ONLY: &
    atmos_im

USE ukca_aero_ctl_mod,             ONLY: &
    drydiam

USE ukca_calc_drydiam_mod,         ONLY: &
    ukca_calc_drydiam

USE ukca_constants,                ONLY: &
    zboltz,                              &
    m_air

USE ukca_radaer_struct_mod,        ONLY: &
    ukca_radaer_struct

USE ukca_volume_mode_mod,          ONLY: &
    ukca_volume_mode

USE um_parvars,                    ONLY: &
    at_extremity

USE um_stashcode_mod,              ONLY: &
    stashcode_glomap_clim_sec

USE yomhook,                       ONLY: &
    lhook,                               &
    dr_hook


IMPLICIT NONE

! Arguments

! Structure for UKCA/radiation interaction
TYPE (ukca_radaer_struct), INTENT(INOUT) :: ukca_radaer

! stashwork54 array
REAL, INTENT(INOUT) :: stashwork54(stash_maxlen(stashcode_glomap_clim_sec,     &
                                                atmos_im))

! Local variables

INTEGER                           :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='PREPARE_FIELDS_FOR_RADAER'

INTEGER :: n_points                   ! Number of points in 3D field
INTEGER :: i                          ! counter
INTEGER :: j                          ! counter
INTEGER :: k                          ! counter
INTEGER :: l                          ! counter
INTEGER :: imode                      ! counter for modes
INTEGER :: icp                        ! counter for components
INTEGER :: section                    ! stash section
INTEGER :: item                       ! stash item
INTEGER :: icode                      ! error code
INTEGER :: ireturn                    ! error code from set_field

INTEGER, PARAMETER   :: im_index = 1  ! internal model index

REAL                 :: mdtmin              ! Minimum value for mdt
REAL, ALLOCATABLE    :: outfield(:,:,:)     ! Array to hold fields for RADAER

LOGICAL, SAVE        :: firstcall = .TRUE.  ! first call to this routine t/f
LOGICAL, ALLOCATABLE :: mask_ndgte(:)       ! Mask for ND threshold

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Find the information for items required from D1 and to copy to stash
IF ( firstcall ) THEN
  CALL glomap_clim_identify_fields()
  firstcall = .FALSE.
END IF

! This has to be appropriate to 85 levels as nmr3d, mmr3d are allocated as 
!  1:model_levels
n_points = row_length*rows*model_levels

IF (.NOT. ALLOCATED(mask_ndgte)) ALLOCATE(mask_ndgte(n_points))
IF (.NOT. ALLOCATED(nd))         ALLOCATE(nd(n_points,nmodes))
IF (.NOT. ALLOCATED(md))         ALLOCATE(md(n_points,nmodes,ncp))
IF (.NOT. ALLOCATED(mdt))        ALLOCATE(mdt(n_points,nmodes))
IF (.NOT. ALLOCATED(aird))       ALLOCATE(aird(n_points))
IF (.NOT. ALLOCATED(pmid_1d))    ALLOCATE(pmid_1d(n_points))
IF (.NOT. ALLOCATED(temp_1d))    ALLOCATE(temp_1d(n_points))
IF (.NOT. ALLOCATED(nmr3d))                                                    &
           ALLOCATE(nmr3d(row_length,rows,model_levels,nmodes))
IF (.NOT. ALLOCATED(mmr3d))                                                    &
           ALLOCATE(mmr3d(row_length,rows,model_levels,nmodes,ncp))

! Copy required fields from D1 array pointers into named item
DO l=1,n_stored_items
  IF (aerofields(l)%required) THEN
    CALL get_d1_fields(l)
  END IF
END DO

IF(.NOT. ALLOCATED( t_theta_levels ) )                                         &
          ALLOCATE( t_theta_levels  ( row_length, rows, model_levels ) )
IF(.NOT. ALLOCATED( pr_theta_levels ) )                                        &
          ALLOCATE( pr_theta_levels ( row_length, rows, model_levels ) )
IF(.NOT. ALLOCATED( qsatmr ) )                                                 &
          ALLOCATE( qsatmr          ( row_length, rows ) )
IF(.NOT. ALLOCATED( qsatmr_wat ) )                                             &
          ALLOCATE( qsatmr_wat      ( row_length, rows ) )

! Prepare 1D arrays for input to lsp_qclear
IF(.NOT. ALLOCATED(q_1d))          ALLOCATE(q_1d(row_length*rows))
IF(.NOT. ALLOCATED(q_clr_1d))      ALLOCATE(q_clr_1d(row_length*rows))
IF(.NOT. ALLOCATED(qsatmr_1d))     ALLOCATE(qsatmr_1d(row_length*rows))
IF(.NOT. ALLOCATED(qsatmr_wat_1d)) ALLOCATE(qsatmr_wat_1d(row_length*rows))
IF(.NOT. ALLOCATED(qcf_1d))        ALLOCATE(qcf_1d(row_length*rows))
IF(.NOT. ALLOCATED(rhcrit_1d))     ALLOCATE(rhcrit_1d(row_length*rows))
IF(.NOT. ALLOCATED(rh_clr_1d))     ALLOCATE(rh_clr_1d(row_length*rows))
IF(.NOT. ALLOCATED(cloud_liq_frac_1d))                                         &
          ALLOCATE(cloud_liq_frac_1d(row_length*rows))
IF(.NOT. ALLOCATED(cloud_blk_frac_1d))                                         &
          ALLOCATE(cloud_blk_frac_1d(row_length*rows))
IF(.NOT. ALLOCATED(rel_humid_frac_clr))                                        &
          ALLOCATE(rel_humid_frac_clr(row_length,rows,model_levels))

! Derive temperature on theta levels
t_theta_levels(:,:,:)= exner_theta_levels(1:row_length,1:rows,1:model_levels)* &
                                    theta(1:row_length,1:rows,1:model_levels)

! Presure on theta levels without halos
pr_theta_levels(:,:,:)  =  p_theta_levels(1:row_length,1:rows,1:model_levels)

DO k=1,model_levels
  CALL qsat_wat_mix ( qsatmr_wat(:,:), t_theta_levels(:,:,k),                  &
                      pr_theta_levels(:,:,k), row_length, rows )
  
  CALL qsat ( qsatmr(:,:), t_theta_levels(:,:,k),                              &
              pr_theta_levels(:,:,k), row_length, rows )

  q_1d              = RESHAPE(q(1:row_length,1:rows,k),                        &
                      (/ SIZE(q(1:row_length,1:rows,k)) /))
  
  qsatmr_1d         = RESHAPE(qsatmr(1:row_length,1:rows),                     &
                      (/ SIZE(qsatmr(1:row_length,1:rows)) /))
  
  qsatmr_wat_1d     = RESHAPE(qsatmr_wat(1:row_length,1:rows),                 &
                      (/ SIZE(qsatmr_wat(1:row_length,1:rows)) /))
  
  qcf_1d            = RESHAPE(qcf(1:row_length,1:rows,k),                      &
                      (/ SIZE(qcf(1:row_length,1:rows,k)) /))
  
  cloud_liq_frac_1d = RESHAPE(cf_liquid(1:row_length,1:rows,k),                &
                      (/ SIZE(cf_liquid(1:row_length,1:rows,k)) /))
  
  cloud_blk_frac_1d = RESHAPE(cf_bulk(1:row_length,1:rows,k),                  &
                      (/ SIZE(cf_bulk(1:row_length,1:rows,k)) /))
  
  rhcrit_1d(:)      = rhcrit(k)
  
  ! Calculate clear sky relative humidity
  CALL lsp_qclear(q_1d, qsatmr_1d, qsatmr_wat_1d,                              &
                  qcf_1d, cloud_liq_frac_1d, cloud_blk_frac_1d,                &
                  rhcrit_1d, q_clr_1d, SIZE(q_1d))
  
  rh_clr_1d(:) = q_clr_1d(:) / qsatmr_wat_1d(:)
  
  rel_humid_frac_clr(1:row_length,1:rows,k) =                                  &
                                         RESHAPE(rh_clr_1d,(/row_length,rows/))
END DO

IF (ALLOCATED(q_1d)) DEALLOCATE(q_1d)
IF (ALLOCATED(q_clr_1d)) DEALLOCATE(q_clr_1d)
IF (ALLOCATED(qsatmr_1d)) DEALLOCATE(qsatmr_1d)
IF (ALLOCATED(qsatmr_wat_1d)) DEALLOCATE(qsatmr_wat_1d)
IF (ALLOCATED(qcf_1d)) DEALLOCATE(qcf_1d)
IF (ALLOCATED(cloud_liq_frac_1d)) DEALLOCATE(cloud_liq_frac_1d)
IF (ALLOCATED(cloud_blk_frac_1d)) DEALLOCATE(cloud_blk_frac_1d)
IF (ALLOCATED(rhcrit_1d)) DEALLOCATE(rhcrit_1d)
IF (ALLOCATED(rh_clr_1d)) DEALLOCATE(rh_clr_1d)
IF (ALLOCATED(qsatmr)) DEALLOCATE(qsatmr)
IF (ALLOCATED(qsatmr_wat)) DEALLOCATE(qsatmr_wat)

! Convert these to 1d for glomap_volume_mode
IF (.NOT. ALLOCATED(rh_clr_1d)) ALLOCATE(rh_clr_1d(n_points))
IF (.NOT. ALLOCATED(q_1d)) ALLOCATE(q_1d(n_points))

rh_clr_1d = RESHAPE( rel_humid_frac_clr(1:row_length, 1:rows, 1:model_levels), &
                  (/ n_points /) )
temp_1d   = RESHAPE( t_theta_levels(1:row_length, 1:rows, 1:model_levels),     &
                  (/ n_points /) )
q_1d      = RESHAPE( q(1:row_length, 1:rows, 1:model_levels), (/ n_points /) )

WHERE (rh_clr_1d < 0.0)
  rh_clr_1d=0.0           ! remove negatives
END WHERE
WHERE (rh_clr_1d > 0.999)
  rh_clr_1d=0.999         ! remove values >= 1.0
END WHERE

pmid_1d(:) = RESHAPE(pr_theta_levels(1:row_length,1:rows,1:model_levels),      &
                                                                 (/n_points/) )

! air density
aird(:) = pmid_1d(:)/(temp_1d(:)*(zboltz*1000000.0))  ! no conc of air (/cm3)

! Allocate arrays for fields to be sent to RADAER

IF (.NOT. ALLOCATED(drydp))    ALLOCATE(drydp(n_points,nmodes))
IF (.NOT. ALLOCATED(wetdp))    ALLOCATE(wetdp(n_points,nmodes))
IF (.NOT. ALLOCATED(rhopar))   ALLOCATE(rhopar(n_points,nmodes))
IF (.NOT. ALLOCATED(dvol))     ALLOCATE(dvol(n_points,nmodes))
IF (.NOT. ALLOCATED(wvol))     ALLOCATE(wvol(n_points,nmodes))
IF (.NOT. ALLOCATED(mdwat))    ALLOCATE(mdwat(n_points,nmodes))
IF (.NOT. ALLOCATED(pvol_wat)) ALLOCATE(pvol_wat(n_points,nmodes))
IF (.NOT. ALLOCATED(pvol))     ALLOCATE(pvol(n_points,nmodes,ncp))

! Fill the 1-D md and nd arrays
nd(:,:) = 0.0
md(:,:,:) = 0.0
mask_ndgte = .FALSE.
DO imode=1,nmodes
  IF (mode(imode)) THEN
    WHERE(nmr3d(:,:,:,imode) < 0.0)
      nmr3d(:,:,:,imode) = 0.0
    END WHERE
    nd(:,imode) = RESHAPE(nmr3d(1:row_length,1:rows,1:model_levels,imode),     &
                          (/ n_points /) )
    nd(:,imode) = nd(:,imode)*aird(:)
    mask_ndgte(:) = (nd(:,imode) > num_eps(imode))
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        WHERE (mask_ndgte(:))
          md(:,imode,icp) = RESHAPE(mmr3d(1:row_length,1:rows,1:model_levels,  &
                                          imode,icp), (/ n_points /) )
          md(:,imode,icp) = md(:,imode,icp)*(m_air/mm(icp))*aird(:)/nd(:,imode)
        ELSEWHERE
          md(:,imode,icp) = mmid(imode)*mfrac_0(imode,icp)
        END WHERE
      END IF
    END DO
  END IF
END DO

! Set total mass array MDT from SUM over individual component MDs
mdt(:,:) = 0.0
DO imode=1,nmodes
  IF (mode(imode)) THEN
    DO icp=1,ncp
      IF (component(imode,icp)) THEN
        mdt(:,imode) = mdt(:,imode) + md(:,imode,icp)
      END IF
    END DO
  END IF
END DO

! Set ND -> 0 where MDT too low and set MDT -> MMID
DO imode=1,nmodes
  IF (mode(imode)) THEN
    mdtmin = mlo(imode)*0.001   ! set equiv. to DPLIM0*0.1
    WHERE (mdt(:,imode) < mdtmin)
      nd(:,imode) = 0.0
      mdt(:,imode) = mmid(imode)
    END WHERE
  END IF
END DO

! Calculate the dry diameters and volumes
CALL ukca_calc_drydiam(n_points, nd, md, mdt, drydp, dvol)

IF ( l_glomap_clim_arg_act ) THEN
  IF ( .NOT. ALLOCATED(drydiam) )                                              &
                  ALLOCATE( drydiam( row_length, rows, model_levels, nmodes ) )
  drydiam(:,:,:,:) = 0.0
  DO imode=1,nmodes
    IF (mode(imode)) THEN
      drydiam(:,:,:,imode) = RESHAPE(drydp(:,imode),                           &
                                             (/row_length,rows,model_levels/) )
    END IF
  END DO
END IF

! Calculate wet diameters, densities, partial volumes...
CALL ukca_volume_mode( n_points, nd, md, mdt, rh_clr_1d, wvol, wetdp,          &
                       rhopar, dvol, drydp, mdwat, pvol,                       &
                       pvol_wat, temp_1d, pmid_1d, q_1d)

! Create STASHwork entries for the section 54 fields required by RADAER

IF (.NOT. ALLOCATED(outfield)) ALLOCATE(outfield(row_length,rows,model_levels))

DO l=1,n_stored_items
  IF ( ( aerofields(l)%item /= imdi ) .AND. ( aerofields(l)%put_stash ) .AND.  &
       ( aerofields(l)%section == stashcode_glomap_clim_sec ) ) THEN
    
    item = aerofields(l)%item
    section = aerofields(l)%section
    
    ! Write fields to ukca_radaer structure
    CALL set_field(item, section, row_length, rows, model_levels,              &
                   ukca_radaer, outfield, ireturn)
    
    IF (ireturn /= 0) THEN
      WRITE(umMessage,'(A,I4,A,I4)') 'item: ', item, ' section: ', section
      CALL umPrint(umMessage, src=RoutineName)
      
      errcode = item
      cmessage = 'outfield contains negative values'
      CALL ereport(Modulename//':'//RoutineName,errcode,cmessage)
    END IF
    
    IF (PrintStatus >= PrStatus_Diag) THEN
      WRITE(umMessage,'(A,L1,A,I2,A,I3)') 'Should field be copied to' //       &
                                          ' stashwork array (T/F)? ',          &
                                          sf(item,section),                    &
                                          ' Section: ', section,               &
                                          ' Item: ', item
      CALL umPrint(umMessage, src=RoutineName)
    END IF
    
    ! Copy requested fields into STASHwork array
    IF (sf(item,section)) THEN
      
      icode = 0
      
      ! DEPENDS ON: copydiag_3d
      CALL copydiag_3d(stashwork54(si(item,section,im_index)),                 &
                       outfield(:,:,:),                                        &
                       tdims%i_len, tdims%j_len, model_levels,                 &
                       0,0,0,0, at_extremity,                                  &
                       stlist(1,stindex(1,item,section,im_index)),len_stlist,  &
                       stash_levels,num_stash_levels+1,                        &
                       atmos_im,section,item,icode,cmessage)
      
      IF (icode >  0) THEN
        errcode = section*1000 + item
        
        WRITE(umMessage,'(A,I6)')'Error from copydiag_3d, errcode: ', errcode
        CALL umPrint(umMessage, src=RoutineName)
        
        CALL ereport(Modulename//':'//RoutineName,errcode,cmessage)
      END IF
    
    END IF   ! IF sf(item,section)
  END IF
END DO       ! 1, n_stored_items

IF (ALLOCATED(outfield))            DEALLOCATE(outfield)
IF (ALLOCATED(mask_ndgte))          DEALLOCATE(mask_ndgte)
IF (ALLOCATED(mmr3d))               DEALLOCATE(mmr3d)
IF (ALLOCATED(nmr3d))               DEALLOCATE(nmr3d)
IF (ALLOCATED(nd))                  DEALLOCATE(nd)
IF (ALLOCATED(md))                  DEALLOCATE(md)
IF (ALLOCATED(mdt))                 DEALLOCATE(mdt)
IF (ALLOCATED(drydp))               DEALLOCATE(drydp)
IF (ALLOCATED(wetdp))               DEALLOCATE(wetdp)
IF (ALLOCATED(rhopar))              DEALLOCATE(rhopar)
IF (ALLOCATED(dvol))                DEALLOCATE(dvol)
IF (ALLOCATED(aird))                DEALLOCATE(aird)
IF (ALLOCATED(wvol))                DEALLOCATE(wvol)
IF (ALLOCATED(mdwat))               DEALLOCATE(mdwat)
IF (ALLOCATED(pvol_wat))            DEALLOCATE(pvol_wat)
IF (ALLOCATED(pvol))                DEALLOCATE(pvol)
IF (ALLOCATED(pr_theta_levels))     DEALLOCATE(pr_theta_levels)
IF (ALLOCATED(t_theta_levels))      DEALLOCATE(t_theta_levels)
IF (ALLOCATED(pmid_1d))             DEALLOCATE(pmid_1d)
IF (ALLOCATED(q_1d))                DEALLOCATE(q_1d)
IF (ALLOCATED(rh_clr_1d))           DEALLOCATE(rh_clr_1d)
IF (ALLOCATED(temp_1d))             DEALLOCATE(temp_1d)
IF (ALLOCATED(rel_humid_frac_clr))  DEALLOCATE(rel_humid_frac_clr)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE prepare_fields_for_radaer

! ##############################################################################

SUBROUTINE get_d1_fields(n)

! To fill the required arrays and using the D1 array pointers

USE glomap_clim_identify_fields_mod, ONLY: &
    aerofields

USE nlsizes_namelist_mod,   ONLY: &
    row_length,                   &
    rows,                         &
    model_levels,                 &
    len_tot

USE parkind1,               ONLY: &
    jprb,                         &
    jpim

USE um_stashcode_mod,       ONLY: &
    stashcode_glomap_clim_sec

USE yomhook,                ONLY: &
    lhook,                        &
    dr_hook


IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN) :: n                       ! index to aerofields array

! Local Variables

INTEGER                           :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='GET_D1_FIELDS'

INTEGER, PARAMETER :: first_mode_item = 101 ! First item for mixing ratio
INTEGER, PARAMETER :: last_mode_item  = 150 ! Final item for mixing ratio

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! GLOMAP-mode aerosol climatology fields for number and mass mixing ratios
IF ( aerofields(n)%section == stashcode_glomap_clim_sec .AND.                  &
     aerofields(n)%item    >= first_mode_item           .AND.                  &
     aerofields(n)%item    <= last_mode_item )          THEN
  
  SELECT CASE(aerofields(n)%item)
  CASE(101)
    nmr3d(:,:,:,mode_nuc_sol)         = gc_nd_nuc_sol(:,:,1:model_levels)
  CASE(102)
    mmr3d(:,:,:,mode_nuc_sol,cp_su)   = gc_nuc_sol_su(:,:,1:model_levels)
  CASE(126)
    mmr3d(:,:,:,mode_nuc_sol,cp_oc)   = gc_nuc_sol_oc(:,:,1:model_levels)
  CASE(103)
    nmr3d(:,:,:,mode_ait_sol)         = gc_nd_ait_sol(:,:,1:model_levels)
  CASE(104)
    mmr3d(:,:,:,mode_ait_sol,cp_su)   = gc_ait_sol_su(:,:,1:model_levels)
  CASE(105)
    mmr3d(:,:,:,mode_ait_sol,cp_bc)   = gc_ait_sol_bc(:,:,1:model_levels)
  CASE(106)
    mmr3d(:,:,:,mode_ait_sol,cp_oc)   = gc_ait_sol_oc(:,:,1:model_levels)
  CASE(107)
    nmr3d(:,:,:,mode_acc_sol)         = gc_nd_acc_sol(:,:,1:model_levels)
  CASE(108)
    mmr3d(:,:,:,mode_acc_sol,cp_su)   = gc_acc_sol_su(:,:,1:model_levels)
  CASE(109)
    mmr3d(:,:,:,mode_acc_sol,cp_bc)   = gc_acc_sol_bc(:,:,1:model_levels)
  CASE(110)
    mmr3d(:,:,:,mode_acc_sol,cp_oc)   = gc_acc_sol_oc(:,:,1:model_levels)
  CASE(111)
    mmr3d(:,:,:,mode_acc_sol,cp_cl)   = gc_acc_sol_ss(:,:,1:model_levels)
  ! CASE(112) gc_acc_sol_du does not exist yet
  CASE(113)
    nmr3d(:,:,:,mode_cor_sol)         = gc_nd_cor_sol(:,:,1:model_levels)
  CASE(114)
    mmr3d(:,:,:,mode_cor_sol,cp_su)   = gc_cor_sol_su(:,:,1:model_levels)
  CASE(115)
    mmr3d(:,:,:,mode_cor_sol,cp_bc)   = gc_cor_sol_bc(:,:,1:model_levels)
  CASE(116)
    mmr3d(:,:,:,mode_cor_sol,cp_oc)   = gc_cor_sol_oc(:,:,1:model_levels)
  CASE(117)
    mmr3d(:,:,:,mode_cor_sol,cp_cl)   = gc_cor_sol_ss(:,:,1:model_levels)
  ! CASE(118) gc_cor_sol_du does not exist yet
  CASE(119)
    nmr3d(:,:,:,mode_ait_insol)       = gc_nd_ait_ins(:,:,1:model_levels)
  CASE(120)
    mmr3d(:,:,:,mode_ait_insol,cp_bc) = gc_ait_ins_bc(:,:,1:model_levels)
  CASE(121)
    mmr3d(:,:,:,mode_ait_insol,cp_oc) = gc_ait_ins_oc(:,:,1:model_levels)
  ! CASE(122) gc_nd_acc_ins does not exist yet
  ! CASE(123) gc_acc_ins_du does not exist yet
  ! CASE(124) gc_nd_cor_ins does not exist yet
  ! CASE(125) gc_cor_ins_du does not exist yet
  CASE DEFAULT
    cmessage=' Item not found in GLOMAP-mode aerosol climatology case statement'
    CALL ereport(RoutineName,aerofields(n)%item,cmessage)
  END SELECT
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE get_d1_fields

! #############################################################################

SUBROUTINE set_field ( item, section, row_length, rows, model_levels,          &
                       ukca_radaer, outfield, ireturn)

! To fill outfield array with correct 3D field using section, item

USE glomap_clim_option_mod,  ONLY: &
    i_glomap_clim_setup,           &
    i_gc_sussocbc_5mode

USE parkind1,                ONLY: &
    jprb,                          &
    jpim

USE ukca_radaer_struct_mod,  ONLY: &
    ukca_radaer_struct

USE um_stashcode_mod,        ONLY: &
    stashcode_glomap_clim_sec,     &
    stashcode_gc_dryd_ait_sol,     &
    stashcode_gc_dryd_acc_sol,     &
    stashcode_gc_dryd_cor_sol,     &
    stashcode_gc_dryd_ait_ins,     &
    stashcode_gc_dryd_acc_ins,     &
    stashcode_gc_dryd_cor_ins,     &
    stashcode_gc_wetd_ait_sol,     &
    stashcode_gc_wetd_acc_sol,     &
    stashcode_gc_wetd_cor_sol,     &
    stashcode_gc_rho_ait_sol,      &
    stashcode_gc_rho_acc_sol,      &
    stashcode_gc_rho_cor_sol,      &
    stashcode_gc_rho_ait_ins,      &
    stashcode_gc_rho_acc_ins,      &
    stashcode_gc_rho_cor_ins,      &
    stashcode_gc_pvol_ait_su_sol,  &
    stashcode_gc_pvol_ait_bc_sol,  &
    stashcode_gc_pvol_ait_oc_sol,  &
    stashcode_gc_pvol_ait_so_sol,  &
    stashcode_gc_pvol_ait_h2o_sol, &
    stashcode_gc_pvol_acc_su_sol,  &
    stashcode_gc_pvol_acc_bc_sol,  &
    stashcode_gc_pvol_acc_oc_sol,  &
    stashcode_gc_pvol_acc_ss_sol,  &
    stashcode_gc_pvol_acc_du_sol,  &
    stashcode_gc_pvol_acc_so_sol,  &
    stashcode_gc_pvol_acc_h2o_sol, &
    stashcode_gc_pvol_cor_su_sol,  &
    stashcode_gc_pvol_cor_bc_sol,  &
    stashcode_gc_pvol_cor_oc_sol,  &
    stashcode_gc_pvol_cor_ss_sol,  &
    stashcode_gc_pvol_cor_du_sol,  &
    stashcode_gc_pvol_cor_so_sol,  &
    stashcode_gc_pvol_cor_h2o_sol, &
    stashcode_gc_pvol_ait_bc_ins,  &
    stashcode_gc_pvol_ait_oc_ins,  &
    stashcode_gc_pvol_acc_du_ins,  &
    stashcode_gc_pvol_cor_du_ins

USE yomhook,                 ONLY: &
    lhook,                         &
    dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN)                      :: item               ! Stash item
INTEGER, INTENT(IN)                      :: section            ! Stash section
INTEGER, INTENT(IN)                      :: row_length         ! field dimension
INTEGER, INTENT(IN)                      :: rows               ! field dimension
INTEGER, INTENT(IN)                      :: model_levels       ! field dimension
TYPE (ukca_radaer_struct), INTENT(INOUT) :: ukca_radaer
REAL,    INTENT(OUT) :: outfield(row_length,rows,model_levels) ! output field
INTEGER, INTENT(OUT) :: ireturn                                ! success/fail

INTEGER, PARAMETER                :: msect = stashcode_glomap_clim_sec * 1000

INTEGER                           :: errcode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='SET_FIELD'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

outfield = -0.1
ireturn = 0
errcode = 0

IF (section == stashcode_glomap_clim_sec) THEN
  SELECT CASE (item)
! Dry diameter
  CASE (stashcode_gc_dryd_ait_sol - msect)        ! ait_sol
    outfield = RESHAPE(drydp(:,2),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%dry_diam(:,:,:,1)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_dryd_acc_sol - msect)        ! acc_sol
    outfield = RESHAPE(drydp(:,3),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%dry_diam(:,:,:,2)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_dryd_cor_sol - msect)        ! cor_sol
    outfield = RESHAPE(drydp(:,4),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%dry_diam(:,:,:,3)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_dryd_ait_ins - msect)        ! ait_ins
    outfield = RESHAPE(drydp(:,5),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%dry_diam(:,:,:,4)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_dryd_acc_ins - msect)        ! acc_ins
    outfield = RESHAPE(drydp(:,6),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_dryd_cor_ins - msect)        ! cor_ins
    outfield = RESHAPE(drydp(:,7),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
! Wet diameter
  CASE (stashcode_gc_wetd_ait_sol - msect)        ! ait_sol
    outfield = RESHAPE(wetdp(:,2),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%wet_diam(:,:,:,1)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_wetd_acc_sol - msect)        ! acc_sol
    outfield = RESHAPE(wetdp(:,3),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%wet_diam(:,:,:,2)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_wetd_cor_sol - msect)        ! cor_sol
    outfield = RESHAPE(wetdp(:,4),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%wet_diam(:,:,:,3)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
! Aerosol density
  CASE (stashcode_gc_rho_ait_sol - msect)         ! ait_sol
    outfield = RESHAPE(rhopar(:,2),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%modal_rho(:,:,:,1)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_rho_acc_sol - msect)         ! acc_sol
    outfield = RESHAPE(rhopar(:,3),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%modal_rho(:,:,:,2)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_rho_cor_sol - msect)         ! cor_sol
    outfield = RESHAPE(rhopar(:,4),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%modal_rho(:,:,:,3)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_rho_ait_ins - msect)         ! ait_ins
    outfield = RESHAPE(rhopar(:,5),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%modal_rho(:,:,:,4)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_rho_acc_ins - msect)         ! acc_ins
    outfield = RESHAPE(rhopar(:,6),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_rho_cor_ins - msect)         ! cor_ins
    outfield = RESHAPE(rhopar(:,7),(/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
! Fractional volume of water in each mode
  CASE (stashcode_gc_pvol_ait_h2o_sol - msect)    ! ait_sol, H2O
    outfield = RESHAPE(pvol_wat(:,mode_ait_sol),                               &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%modal_wtv(:,:,:,1)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_acc_h2o_sol - msect)    ! acc_sol, H2O
    outfield = RESHAPE(pvol_wat(:,mode_acc_sol),                               &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%modal_wtv(:,:,:,2)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_cor_h2o_sol - msect)    ! cor_sol, H2O
    outfield = RESHAPE(pvol_wat(:,mode_cor_sol),                               &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%modal_wtv(:,:,:,3)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
! Partial volume ait_sol
  CASE (stashcode_gc_pvol_ait_su_sol - msect)     ! ait_sol, SU
    outfield = RESHAPE(pvol(:,mode_ait_sol,cp_su),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,1)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_ait_bc_sol - msect)     ! ait_sol, BC
    outfield = RESHAPE(pvol(:,mode_ait_sol,cp_bc),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,2)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_ait_oc_sol - msect)     ! ait_sol, OC
    outfield = RESHAPE(pvol(:,mode_ait_sol,cp_oc),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,3)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
                                                  ! ait_sol, SS does not exist
                                                  ! ait_sol, DU does not exist
  CASE (stashcode_gc_pvol_ait_so_sol - msect)     ! ait_sol, SO
    outfield = RESHAPE(pvol(:,mode_ait_sol,cp_so),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
! Partial volume acc_sol
  CASE (stashcode_gc_pvol_acc_su_sol - msect)     ! acc_sol, SU
    outfield = RESHAPE(pvol(:,mode_acc_sol,cp_su),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,4)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_acc_bc_sol - msect)     ! acc_sol, BC
    outfield = RESHAPE(pvol(:,mode_acc_sol,cp_bc),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,5)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_acc_oc_sol - msect)     ! acc_sol, OC
    outfield = RESHAPE(pvol(:,mode_acc_sol,cp_oc),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,6)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_acc_ss_sol - msect)     ! acc_sol, SS
    outfield = RESHAPE(pvol(:,mode_acc_sol,cp_cl),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,7)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_acc_du_sol - msect)     ! acc_sol, DU
    outfield = RESHAPE(pvol(:,mode_acc_sol,cp_du),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_acc_so_sol - msect)     ! acc_sol, SO
    outfield = RESHAPE(pvol(:,mode_acc_sol,cp_so),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
! Partial volume cor_sol
  CASE (stashcode_gc_pvol_cor_su_sol - msect)     ! cor_sol, SU
    outfield = RESHAPE(pvol(:,mode_cor_sol,cp_su),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,8)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_cor_bc_sol - msect)     ! cor_sol, BC
    outfield = RESHAPE(pvol(:,mode_cor_sol,cp_bc),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,9)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_cor_oc_sol - msect)     ! cor_sol, OC
    outfield = RESHAPE(pvol(:,mode_cor_sol,cp_oc),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,10) = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_cor_ss_sol - msect)     ! cor_sol, SS
    outfield = RESHAPE(pvol(:,mode_cor_sol,cp_cl),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,11)  = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_cor_du_sol - msect)     ! cor_sol, DU
    outfield = RESHAPE(pvol(:,mode_cor_sol,cp_du),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_cor_so_sol - msect)     ! cor_sol, SO
    outfield = RESHAPE(pvol(:,mode_cor_sol,cp_so),                             &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
! Partial volume ait_ins
                                                  ! ait_ins, SU does not exist
  CASE (stashcode_gc_pvol_ait_bc_ins - msect)     ! ait_ins, BC
    outfield = RESHAPE(pvol(:,mode_ait_insol,cp_bc),                           &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,12) = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
  CASE (stashcode_gc_pvol_ait_oc_ins - msect)     ! ait_ins, OC
    outfield = RESHAPE(pvol(:,mode_ait_insol,cp_oc),                           &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE (i_gc_sussocbc_5mode)
      ukca_radaer%comp_vol(:,:,:,13) = outfield
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
                                                  ! ait_ins, SS does not exist
                                                  ! ait_ins, DU does not exist
                                                  ! ait_ins, SO does not exist
    
! Partial volume acc_ins
   
                                                  ! acc_ins, SU does not exist
                                                  ! acc_ins, BC does not exist
                                                  ! acc_ins, OC does not exist
                                                  ! acc_ins, SS does not exist
  
  CASE (stashcode_gc_pvol_acc_du_ins - msect)     ! acc_ins, DU
    outfield = RESHAPE(pvol(:,mode_acc_insol,cp_du),                           &
                       (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
                                                  ! acc_ins, SO does not exist
    
! Partial volume cor_ins
                                                  ! cor_ins, SU does not exist
                                                  ! cor_ins, BC does not exist
                                                  ! cor_ins, OC does not exist
                                                  ! cor_ins, SS does not exist
  
  CASE (stashcode_gc_pvol_cor_du_ins - msect)     ! cor_ins, DU
     outfield = RESHAPE(pvol(:,mode_cor_insol,cp_du),                          &
                        (/row_length,rows,model_levels/))
    
    SELECT CASE(i_glomap_clim_setup)
    CASE DEFAULT
      errcode = ABS(item)
    END SELECT
    
                                                  ! cor_ins, SO does not exist
  
! Default
  CASE DEFAULT
    cmessage = ' Item not found in case statement'
    errcode = ABS(item)
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END SELECT
  
  IF ( errcode /= 0 ) THEN
    cmessage = 'Item should not be used for choice of i_glomap_clim_setup'
    CALL ereport(RoutineName,errcode,cmessage)
  END IF
  
END IF ! section == stashcode_glomap_clim_sec

! Set flag if field not filled
IF (ANY(outfield < 0.0)) ireturn = 1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_field

END MODULE glomap_clim_fields_mod
