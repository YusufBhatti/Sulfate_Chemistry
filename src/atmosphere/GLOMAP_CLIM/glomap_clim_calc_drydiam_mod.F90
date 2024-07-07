! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Populate drydiam field
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! Procedure:
!   1) CALL glomap_clim_identify_fields to identify item numbers required
!
!   2) CALL get_d1_fields to populate the nmr3d and mmr3d fields
!
!   3) CALL ukca_calc_drydiam routine to calculate the drydiam field required 
!       by ACTIVATE
!
! ---------------------------------------------------------------------
MODULE glomap_clim_calc_drydiam_mod

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='GLOMAP_CLIM_CALC_DRYDIAM_MOD'

CONTAINS

SUBROUTINE glomap_clim_calc_drydiam ( t_theta_levels , pr_theta_levels )

USE glomap_clim_fields_mod,          ONLY: &
    get_d1_fields,                         &
    aird,                                  &
    drydp,                                 &
    dvol,                                  &
    md,                                    &
    mdt,                                   &
    mmr3d,                                 &
    nd,                                    &
    nmr3d,                                 &
    pmid_1d,                               &
    temp_1d

USE glomap_clim_identify_fields_mod, ONLY: &
    glomap_clim_identify_fields,           &
    aerofields,                            &
    n_stored_items

USE nlsizes_namelist_mod,            ONLY: &
    row_length,                            &
    rows,                                  &
    model_levels

USE parkind1,                        ONLY: &
    jpim,                                  &
    jprb

USE ukca_aero_ctl_mod,               ONLY: &
    drydiam

USE ukca_calc_drydiam_mod,           ONLY: &
    ukca_calc_drydiam

USE ukca_constants,                  ONLY: &
    zboltz,                                &
    m_air

USE ukca_mode_setup,                 ONLY: &
    component,                             &
    mfrac_0,                               &
    mlo,                                   &
    mm,                                    &
    mmid,                                  &
    mode,                                  &
    ncp,                                   &
    nmodes,                                &
    num_eps

USE yomhook,                         ONLY: &
    lhook,                                 &
    dr_hook

IMPLICIT NONE

! Arguments

! Temperature
REAL, INTENT(IN) :: t_theta_levels(:,:,:)

! Pressure on theta levels without halos
REAL, INTENT(IN) :: pr_theta_levels(:,:,:)

! Local variables

LOGICAL, SAVE        :: firstcall = .TRUE.  ! first call to this routine t/f

LOGICAL, ALLOCATABLE :: mask_ndgte(:)       ! Mask for ND threshold

INTEGER :: n_points                   ! Number of points in 3D field
INTEGER :: imod                       ! counter for modes
INTEGER :: icp                        ! counter for components
INTEGER :: n                          ! counter

REAL    :: mdtmin                     ! Minimum value for mdt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName='GLOMAP_CLIM_CALC_DRYDIAM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

IF (.NOT. ALLOCATED(drydiam))                                                  &
           ALLOCATE(drydiam(row_length,rows,model_levels,nmodes))

! Find the information for items required from D1
IF ( firstcall ) THEN
  CALL glomap_clim_identify_fields()
  firstcall = .FALSE.
END IF

n_points = row_length * rows * model_levels

IF (.NOT. ALLOCATED(temp_1d)) ALLOCATE(temp_1d(n_points))
temp_1d(:) = RESHAPE( t_theta_levels(1:row_length, 1:rows, 1:model_levels),    &
                                                               (/ n_points /) )

IF (.NOT. ALLOCATED(pmid_1d)) ALLOCATE(pmid_1d(n_points))
pmid_1d(:) = RESHAPE( pr_theta_levels(1:row_length, 1:rows, 1:model_levels),   &
                                                               (/ n_points /) )

IF (.NOT. ALLOCATED(aird)) ALLOCATE(aird(n_points))
aird(:) = pmid_1d(:) / ( temp_1d(:) * ( zboltz * 1000000.0 ) )

IF (.NOT. ALLOCATED(nmr3d))                                                    &
           ALLOCATE(nmr3d(row_length,rows,model_levels,nmodes))

IF (.NOT. ALLOCATED(mmr3d))                                                    &
           ALLOCATE(mmr3d(row_length,rows,model_levels,nmodes,ncp))

! Copy required fields from D1 array pointers into nmr3d & mmr3d
DO n=1,n_stored_items
  IF (aerofields(n)%required) THEN
    CALL get_d1_fields(n)
  END IF
END DO

IF (.NOT. ALLOCATED(mask_ndgte)) ALLOCATE(mask_ndgte(n_points))
IF (.NOT. ALLOCATED(nd))         ALLOCATE(nd(n_points,nmodes))
IF (.NOT. ALLOCATED(md))         ALLOCATE(md(n_points,nmodes,ncp))
! Fill the 1-D md and nd arrays
nd(:,:) = 0.0
md(:,:,:) = 0.0
mask_ndgte = .FALSE.
DO imod=1,nmodes
  IF (mode(imod)) THEN
    WHERE(nmr3d(:,:,:,imod) < 0.0)
      nmr3d(:,:,:,imod) = 0.0
    END WHERE
    
    nd(:,imod) = RESHAPE(nmr3d(1:row_length,1:rows,1:model_levels,imod),       &
                                                               (/ n_points /) )
    nd(:,imod) = nd(:,imod)*aird(:)
    mask_ndgte(:) = (nd(:,imod) > num_eps(imod))
    DO icp=1,ncp
      IF (component(imod,icp)) THEN
        WHERE (mask_ndgte(:))
          md(:,imod,icp) = RESHAPE(mmr3d(1:row_length,1:rows,1:model_levels,   &
                                         imod,icp), (/ n_points /) )
          md(:,imod,icp) = md(:,imod,icp)*(m_air/mm(icp))*aird(:)/nd(:,imod)
        ELSEWHERE
          md(:,imod,icp) = mmid(imod)*mfrac_0(imod,icp)
        END WHERE
      END IF
    END DO
  END IF
END DO

IF (.NOT. ALLOCATED(mdt))        ALLOCATE(mdt(n_points,nmodes))
! Set total mass array MDT from SUM over individual component MDs
mdt(:,:) = 0.0
DO imod=1,nmodes
  IF (mode(imod)) THEN
    DO icp=1,ncp
      IF (component(imod,icp)) THEN
        mdt(:,imod) = mdt(:,imod) + md(:,imod,icp)
      END IF
    END DO
  END IF
END DO

! Set ND -> 0 where MDT too low and set MDT -> MMID
DO imod=1,nmodes
  IF (mode(imod)) THEN
    mdtmin = mlo(imod)*0.001   ! set equiv. to DPLIM0*0.1
    WHERE (mdt(:,imod) < mdtmin)
      nd(:,imod) = 0.0
      mdt(:,imod) = mmid(imod)
    END WHERE
  END IF
END DO

IF (.NOT. ALLOCATED(drydp)) ALLOCATE(drydp(n_points,nmodes))
IF (.NOT. ALLOCATED(dvol))  ALLOCATE(dvol(n_points,nmodes))

! Calculate the dry diameters and volumes
CALL ukca_calc_drydiam(n_points, nd, md, mdt, drydp, dvol)

! Initialise to zero
drydiam(:,:,:,:) = 0.0

DO imod=1,nmodes
  IF (mode(imod)) THEN
    drydiam(:,:,:,imod)=RESHAPE(drydp(:,imod),(/row_length,rows,model_levels/))
  END IF
END DO

IF (ALLOCATED(dvol))                DEALLOCATE(dvol)
IF (ALLOCATED(drydp))               DEALLOCATE(drydp)
IF (ALLOCATED(mdt))                 DEALLOCATE(mdt)
IF (ALLOCATED(md))                  DEALLOCATE(md)
IF (ALLOCATED(nd))                  DEALLOCATE(nd)
IF (ALLOCATED(mask_ndgte))          DEALLOCATE(mask_ndgte)
IF (ALLOCATED(mmr3d))               DEALLOCATE(mmr3d)
IF (ALLOCATED(nmr3d))               DEALLOCATE(nmr3d)
IF (ALLOCATED(aird))                DEALLOCATE(aird)
IF (ALLOCATED(pmid_1d))             DEALLOCATE(pmid_1d)
IF (ALLOCATED(temp_1d))             DEALLOCATE(temp_1d)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE glomap_clim_calc_drydiam

END MODULE glomap_clim_calc_drydiam_mod
