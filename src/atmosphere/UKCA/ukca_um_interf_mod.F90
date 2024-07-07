! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!  To gather required fields from D1, and replace them after chemistry.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  Fortran 2003
!
! ----------------------------------------------------------------------
MODULE ukca_um_interf_mod

USE ukca_d1_defs
USE ukca_tracer_stash,   ONLY: a_max_ukcavars
USE d1_array_mod, ONLY: d1_object_type, d1_section, d1_item, d1_address,    &
                        d1_length, d1_grid_type, d1_no_levels, d1_halo_type,&
                        prognostic, d1, d1_addr, no_obj_d1
USE atm_fields_bounds_mod, ONLY: tdims, tdims_s, udims
USE nlsizes_namelist_mod, ONLY: global_row_length, global_rows,             &
                                land_field, row_length, rows, tr_ukca,      &
                                len_tot, model_levels, n_obj_d1_max,        &
                                bl_levels

USE ukca_ntp_mod,         ONLY: ntp_type, dim_ntp, ntp_init, print_all_ntp, &
                                stash2ntpindex, ntp_dealloc, name2ntpindex
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC ::  putd1flds, getd1flds, ukca_pr_inputs, ukca_um_dealloc,           &
           ukca_um_d1_initialise

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_UM_INTERF_MOD'

! arrays filled from D1, allocated dynamically with correct halo read
! in from addressing array
INTEGER, ALLOCATABLE, PUBLIC :: conv_cloud_base(:,:)
INTEGER, ALLOCATABLE, PUBLIC :: conv_cloud_top(:,:)
INTEGER, ALLOCATABLE, PUBLIC :: kent(:,:)
INTEGER, ALLOCATABLE, PUBLIC :: kent_dsc(:,:)

LOGICAL, ALLOCATABLE, PUBLIC :: land_sea_mask(:,:)

! temporary array for reading NTP data from D1
REAL, ALLOCATABLE, PUBLIC :: ntp_data(:,:,:) 

! Non transported prognostics
TYPE(ntp_type), SAVE, PUBLIC :: all_ntp(dim_ntp)

! Real variables in order of allocation in getd1flds
REAL, ALLOCATABLE, PUBLIC :: emission1(:,:)
REAL, ALLOCATABLE, PUBLIC :: all_emissions(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: dust_flux(:,:,:)     ! dust emissions 
                                                  ! from 6 bin scheme
REAL, ALLOCATABLE, PUBLIC :: theta(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: soil_moisture_layer1(:)
REAL, ALLOCATABLE, PUBLIC :: q(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: qcf(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: conv_cloud_lwp(:,:)
REAL, ALLOCATABLE, PUBLIC :: Tstar(:,:)
REAL, ALLOCATABLE, PUBLIC :: zbl(:,:)      ! BL height
REAL, ALLOCATABLE, PUBLIC :: Rough_length(:,:)
REAL, ALLOCATABLE, PUBLIC :: seaice_frac(:,:)
REAL, ALLOCATABLE, PUBLIC :: dms_land_ems(:,:)  ! DMS emiss if
REAL, ALLOCATABLE, PUBLIC :: um_ozone(:,:,:)   ! i.e. from D1
REAL, ALLOCATABLE, PUBLIC :: chloro_sea(:,:)  ! surface chlorophyll
REAL, ALLOCATABLE, PUBLIC :: so4_aitken(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: so4_accum(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: co2_interactive(:,:,:)  ! CO2 MMR

! Aerosol mmr / numbers from CLASSIC for heterogeneous reactions
! (note that so4_aitken & so4_accum can also be used for heterogeneous
! chemistry but they are declared somewhere else)
REAL, ALLOCATABLE, PUBLIC :: soot_fresh    (:,:,:) ! mmr (kg kg-1)
REAL, ALLOCATABLE, PUBLIC :: soot_aged     (:,:,:)
REAL, ALLOCATABLE, PUBLIC :: ocff_fresh    (:,:,:)
REAL, ALLOCATABLE, PUBLIC :: ocff_aged     (:,:,:)

REAL, ALLOCATABLE, PUBLIC :: dms_sea_conc(:,:)  ! no CLASSIC
REAL, ALLOCATABLE, PUBLIC :: vertvel(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: conv_cloud_amount(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: Frac_types(:,:)
REAL, ALLOCATABLE, PUBLIC :: laift_lp(:,:)    ! LAI
REAL, ALLOCATABLE, PUBLIC :: canhtft_lp(:,:)  ! canopy height
REAL, ALLOCATABLE, PUBLIC :: canwctile_lp(:,:) ! canopy WC
REAL, ALLOCATABLE, PUBLIC :: Tstar_tile(:,:)
REAL, ALLOCATABLE, PUBLIC :: z0tile_lp(:,:)
REAL, ALLOCATABLE, PUBLIC :: Snow_tile(:,:)
REAL, ALLOCATABLE, PUBLIC :: rho_r2(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: qcl(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: exner_rho_levels(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: area_cloud_fraction(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: cloud_frac(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: cloud_liq_frac(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: biogenic      (:,:,:)
REAL, ALLOCATABLE, PUBLIC :: fland(:)      ! land fraction

! surf_albedo is interpolated at every timestep from land_albedo_all
! which is calculated on radiation timesteps using the SW fluxes.
REAL, ALLOCATABLE, PUBLIC :: land_albedo(:,:)

REAL, ALLOCATABLE, PUBLIC :: exner_theta_levels(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: p_rho_levels(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: p_theta_levels(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: pstar(:,:)
REAL, ALLOCATABLE, PUBLIC :: net_surf_SW(:,:)
REAL, ALLOCATABLE, PUBLIC :: tot_surf_SW(:,:)
REAL, ALLOCATABLE, PUBLIC :: sea_salt_film (:,:,:) ! number (m-3) 
REAL, ALLOCATABLE, PUBLIC :: sea_salt_jet  (:,:,:)
REAL, ALLOCATABLE, PUBLIC :: sulphate_od(:,:,:)  ! optical depth
REAL, ALLOCATABLE, PUBLIC :: rhokh_mix(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: dtrdz_charney_grid(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: we_lim(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: t_frac(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: zrzi(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: we_lim_dsc(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: t_frac_dsc(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: zrzi_dsc(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: zhsc(:,:)
REAL, ALLOCATABLE, PUBLIC :: U_scalar_10m(:,:)
REAL, ALLOCATABLE, PUBLIC :: surf_hf(:,:)
REAL, ALLOCATABLE, PUBLIC :: stcon(:,:,:)       ! stomatal cond
REAL, ALLOCATABLE, PUBLIC :: u_s(:,:)     ! surf frict velocity
REAL, ALLOCATABLE, PUBLIC :: bl_tke(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: cloud_liq_water(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: ls_rain3d(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: ls_snow3d(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: ls_ppn_frac(:,:,:)

! Melting rate ice crystals kg/kg/s
REAL, ALLOCATABLE, PUBLIC :: ice_melt(:,:,:)

! Melting rate snow aggs. kg/kg/s
REAL, ALLOCATABLE, PUBLIC :: snow_melt(:,:,:)

! Riming rate of ice crystals kg/kg/s
REAL, ALLOCATABLE, PUBLIC :: rim_cry(:,:,:) 

! Riming rate of ice aggregates kg/kg/s  
REAL, ALLOCATABLE, PUBLIC :: rim_agg(:,:,:) 

! Rain Autoconversion rate kg/kg/s
REAL, ALLOCATABLE, PUBLIC :: autoconv(:,:,:)

! Rain Accretion rate kg/kg/s
REAL, ALLOCATABLE, PUBLIC :: accretion(:,:,:)

REAL, ALLOCATABLE, PUBLIC :: conv_rain3d(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: conv_snow3d(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: ch4_wetl_emiss(:,:)
REAL, ALLOCATABLE, PUBLIC :: pv_on_theta_mlevs(:,:,:)
REAL, ALLOCATABLE, PUBLIC :: tropopause_height(:,:)

! Other variables
REAL, ALLOCATABLE, PUBLIC :: diag_temp(:,:,:)   ! To hold data from D1
REAL, ALLOCATABLE, PUBLIC :: trmol_post_atmstep(:,:,:,:)
REAL, ALLOCATABLE, PUBLIC :: mode_diags(:,:,:,:)
REAL, ALLOCATABLE, PUBLIC :: strat_fluxdiags(:,:,:,:)
REAL, ALLOCATABLE, PUBLIC :: totnodens(:,:,:)  ! density in molecs/m^3
REAL, ALLOCATABLE, PUBLIC :: rel_humid_frac(:,:,:)

! sat vap pressure with respect to liquid water irrespective of temperature
REAL, ALLOCATABLE, PUBLIC :: qsvp(:,:,:)

! sat mixing ratio with respect to liquid water irrespective of temperature
REAL, ALLOCATABLE, PUBLIC :: qsmr(:,:,:)

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE ukca_um_d1_initialise( &
m_atm_modl, klev1, first)

USE umPrintMgr, ONLY: PrintStatus, PrStatus_Oper, PrStatus_Diag, &
                      umPrint, umMessage
USE stash_array_mod, ONLY: totitems, stlist
USE stparam_mod,  ONLY: st_d1pos, st_macrotag
USE errormessagelength_mod, ONLY: errormessagelength
USE missing_data_mod, ONLY: imdi
USE ereport_mod,  ONLY: ereport
IMPLICIT NONE

! Arguments
INTEGER :: m_atm_modl
INTEGER :: klev1
LOGICAL :: first

! Internal variables
INTEGER    :: get_fld_type      ! UM function

! d1 settings
INTEGER :: ptd1    
INTEGER :: section 
INTEGER :: item    
INTEGER :: levs    
INTEGER :: length  
INTEGER :: addr    
INTEGER :: halo_typ
INTEGER :: grid_typ
INTEGER :: stashcode
INTEGER :: tag
INTEGER :: field_typ         ! Field type
INTEGER :: i,j, l  ! Loop counters
INTEGER :: errcode ! Error code
CHARACTER(LEN=errormessagelength) :: cmessage

IF (first) THEN
  IF (PrintStatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'(A)') 'UKCA Prognostics and Diagnostics from D1:'
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,A)') 'section,item,levels,length,address,',&
                             'halo_type,grid_type,field_type'
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  END IF
  DO i=1,no_obj_D1(m_atm_modl)
    section = d1_addr(D1_section,i,m_atm_modl)
    item = d1_addr(D1_item,i,m_atm_modl)
    levs = d1_addr(d1_no_levels,i,m_atm_modl)
    length = d1_addr(d1_length,i,m_atm_modl)
    addr    = d1_addr(d1_address,i,m_atm_modl)
    halo_typ = d1_addr(d1_halo_type,i,m_atm_modl)
    grid_typ = d1_addr(d1_grid_type,i,m_atm_modl)
    ! DEPENDS ON: get_fld_type
    field_typ = get_fld_type(grid_typ)
    DO j=1,Nukca_D1items
      IF (UkcaD1Codes(j)%section == section .AND.                &
          UkcaD1Codes(j)%item == item .AND.                      &
          d1_addr(d1_object_type,i,m_atm_modl) == prognostic     &
          .AND. UkcaD1codes(j)%Prognostic) THEN
        UkcaD1Codes(j)%n_levels=levs
        UkcaD1Codes(j)%address=addr
        UkcaD1Codes(j)%length=length
        UkcaD1Codes(j)%halo_type=halo_typ
        UkcaD1Codes(j)%grid_type=grid_typ
        UkcaD1Codes(j)%field_type=field_typ
        IF (PrintStatus >= PrStatus_Diag) THEN
          WRITE(umMessage,'(A,4I5,2I10,3I4)') 'P',j, section,item,&
                          levs,length,addr,halo_typ,grid_typ,field_typ
          CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
        END IF
      END IF
    END DO
  END DO

  !       Diagnostics loop through all stashlist items, check
  !       for tag value, read adress etc into UkcaD1Codes.

  DO l=1,totitems
    tag=stlist(st_macrotag,l)-(1000*(stlist(st_macrotag,l)/1000))
    IF (tag == 98) THEN
      ptd1      = stlist(st_D1pos,l)
      section   = d1_addr(d1_section,  ptd1,m_atm_modl)
      item      = d1_addr(d1_item,     ptd1,m_atm_modl)
      levs      = d1_addr(d1_no_levels,ptd1,m_atm_modl)
      length    = d1_addr(d1_length,   ptd1,m_atm_modl)
      addr      = d1_addr(d1_address,  ptd1,m_atm_modl)
      halo_typ  = d1_addr(d1_halo_type,ptd1,m_atm_modl)
      grid_typ  = d1_addr(d1_grid_type,ptd1,m_atm_modl)
      ! DEPENDS ON: get_fld_type
      field_typ = get_fld_type(grid_typ)
      stashcode = section*1000 + item
      DO j=1,Nukca_D1items
        IF (UkcaD1Codes(j)%section == section .AND.              &
            UkcaD1Codes(j)%item    == item    .AND.              &
            .NOT. UkcaD1Codes(j)%Prognostic) THEN
          UkcaD1Codes(j)%n_levels  = levs
          UkcaD1Codes(j)%address   = addr
          UkcaD1Codes(j)%length    = length
          UkcaD1Codes(j)%halo_type = halo_typ
          UkcaD1Codes(j)%grid_type = grid_typ
          UkcaD1Codes(j)%field_type= field_typ
          IF (PrintStatus >= PrStatus_Diag) THEN
            WRITE(umMessage,'(A,4I5,I10,3I4)') 'D',section,item, &
                         levs,length,addr,halo_typ,grid_typ,field_typ
            CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
          END IF
        END IF
      END DO
    END IF  ! tag
  END DO    ! L

END IF ! first call?
! ----------------------------------------------------------------------
! Get data from D1
! ----------------------------------------------------------------------
! Check if all items selected have been identified in D1.
! If any are missing abort. If all selected items have been 
! identified then call getd1flds for all items required.
! Two seperate loops allow all missing inputs to be indentified
! without having to run the code multiple times.
! ----------------------------------------------------------------------
! Do checking
errcode = 0
DO i=1,Nukca_D1items
  IF (UkcaD1codes(i)%address == imdi .AND.                       &
      UkcaD1codes(i)%required) THEN
    errcode  = errcode + 1
    cmessage = 'Item address not found in D1 array:  '
    WRITE(umMessage,'(A37,I5)') cmessage (1:37),                         &
      UkcaD1Codes(i)%section*1000 + UkcaD1Codes(i)%item
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  END IF
END DO
IF (errcode > 0) THEN
  cmessage = 'Some item addresses not found in D1 array'
  CALL ereport('UKCA_MAIN', errcode, cmessage)
END IF

! Copy fields from D1 array into named item, allocate arrays

DO i=1,Nukca_D1items
  IF (UkcaD1Codes(i)%required) THEN
    CALL getd1flds( &
i, klev1, first)
  END IF
END DO


END SUBROUTINE ukca_um_d1_initialise
! ----------------------------------------------------------------------
! Get prognostic and diagnostic fields from D1
! ----------------------------------------------------------------------
SUBROUTINE getd1flds( &
n, klev1, first)

USE set_rad_steps_mod,         ONLY: l_rad_step_prog
USE ukca_extract_d1_data_mod,  ONLY: ukca_extract_d1_data
USE ukca_set_array_bounds_mod, ONLY: ukca_set_array_bounds
USE ukca_tracer_stash,         ONLY: a_max_ukcavars
USE errormessagelength_mod,    ONLY: errormessagelength
USE ereport_mod,               ONLY: ereport
USE umPrintMgr

USE dust_parameters_mod,    ONLY: ndiv

IMPLICIT NONE

INTEGER, INTENT(IN)     :: n             ! id of array
INTEGER, INTENT(IN)     :: klev1         ! First Level for T/Tr/Q grid
LOGICAL, INTENT(IN)     :: first         ! first call?

INTEGER, SAVE           :: next_tr       ! count tracers
INTEGER, SAVE           :: idust         ! count dust divs
INTEGER, SAVE           :: next_em       ! count emissions
INTEGER                 :: i1,i2,j1,j2   ! array bounds
INTEGER :: varindex

REAL, ALLOCATABLE       :: rtemp_2d(:,:) ! Dummy array

! ErrorStatus
INTEGER                    :: errcode=0     ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)         :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GETD1FLDS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The index n is used to look up the STASH section and item for the
! data which is being retrieved from D1. 
! The subroutine gets the appropriate data for this section/item
! from the D1 array and allocates an array to hold the data

! Attempt to match the section and item requested in turn, with one 
! large IF block

! Start with prognostic variables - all in section ukca_sect
IF (UkcaD1codes(n)%section == UKCA_sect) THEN

! If in UKCA section and not a tracer, then it is a non-transported prognostic
! Non transported prognostic variables - each have their own 
! 3D array held in a single vector of ntp_type. 

    ! For ENDGame compatibility - extract all levels (klev1=0 or 1) 
    ! but use only 1:model_levels 

    ! First get data to a temporary 3D array ntp_data with all levels
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    IF (.NOT. ALLOCATED(ntp_data)) THEN
      ALLOCATE(ntp_data(i1:i2,j1:j2,                           &
               klev1:ukcaD1codes(n)%len_dim3))
    END IF
    CALL ukca_extract_d1_data(                                   &
     first,n,ntp_data)

    ! Now pass this into the all_ntp structure. Don't include level 0
    varindex = stash2ntpindex(all_ntp, UkcaD1codes(n)%section,                 &
      UkcaD1codes(n)%item)
    ALLOCATE(all_ntp(varindex)%data_3d(i1:i2,j1:j2,ukcaD1codes(n)%len_dim3))
    all_ntp(varindex)%data_3d(:,:,1:model_levels) =                            &
      ntp_data(:,:,1:model_levels)

  ! Using new emission system - set up dummy arrays
  IF (.NOT. ALLOCATED(emission1)) THEN
    ALLOCATE (emission1     (1,1))
    ALLOCATE (all_emissions (1,1,1))
  END IF
  IF (.NOT. ALLOCATED(em_index))                               &
    ALLOCATE (em_index(1))
  emission1     (:,:)   = 0.0
  all_emissions (:,:,:) = 0.0
  em_index      (:)     = 0

  ! Add  Dust emissions
ELSE IF (UkcaD1codes(n)%section == 3 .AND.                       &
        UkcaD1codes(n)%item >= 401 .AND.                         &
        UkcaD1codes(n)%item <= 401+ndiv-1 ) THEN
  IF (.NOT. ALLOCATED(dust_flux)) THEN
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(dust_flux(i1:i2,j1:j2,ndiv))
    idust = 1
  END IF
  CALL ukca_extract_d1_data(                                     &
     first,n,dust_flux(:,:,idust))
  idust = idust + 1

  !     Prognostics: fill appropriate array

ELSE IF (UkcaD1codes(n)%section == 0 .AND.                       &
        UkcaD1codes(n)%prognostic) THEN
  SELECT CASE(UkcaD1codes(n)%item)
  CASE (4)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(theta(i1:i2,j1:j2,                                  &
                        klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,theta)
  CASE (9)
! Soil moisture - Read data for all soil levels, but use only level 1
    ALLOCATE(rtemp_2d(ukcaD1codes(n)%len_dim1,ukcaD1codes(n)%len_dim2))
    ALLOCATE( soil_moisture_layer1(ukcaD1codes(n)%len_dim1) )
    CALL ukca_extract_d1_data(                                   &
     first,n,rtemp_2d)
    soil_moisture_layer1(:) = rtemp_2d(:,1)
    DEALLOCATE(rtemp_2d) 
  CASE (10)
     ! if Q has been allocated then this does not need to be taken again
     ! and will in fact cause an error if it is.
    IF (.NOT. ALLOCATED(q)) THEN
      CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
      ALLOCATE(q(i1:i2,j1:j2,                                   &
           klev1:ukcaD1codes(n)%len_dim3))
      CALL ukca_extract_d1_data(                                   &
       first,n,q)

    END IF
  CASE (12)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(qcf(i1:i2,j1:j2,                                    &
                        klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,qcf)
  CASE (16)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(conv_cloud_lwp(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,conv_cloud_lwp)
  CASE (24)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(tstar(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,tstar)
  CASE (25)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(zbl(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,zbl)
  CASE (26)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(Rough_length(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,Rough_length)
  CASE (30)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(land_sea_mask(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,land_sea_mask)
  CASE (31)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(seaice_frac(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,seaice_frac)
  CASE (59)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(dms_land_ems(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,dms_land_ems)
  CASE (60)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(um_ozone(i1:i2,j1:j2,                               &
                       klev1:ukcaD1codes(n)%len_dim3))
    um_ozone = RESHAPE(d1(UkcaD1Codes(n)%address:                &
       UkcaD1Codes(n)%address+UkcaD1Codes(n)%length-1),          &
       (/SIZE(um_ozone,dim=1),SIZE(um_ozone,dim=2),              &
       SIZE(um_ozone,dim=3)/))
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
  CASE (96)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(chloro_sea(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,chloro_sea)
  CASE (103)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(so4_aitken(i1:i2,j1:j2,                             &
                        klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,so4_aitken)
  CASE (104)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(so4_accum(i1:i2,j1:j2,                              &
                        klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,so4_accum)
  CASE (108)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(soot_fresh(i1:i2,j1:j2,                             &
                        klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,soot_fresh)
  CASE (109)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(soot_aged(i1:i2,j1:j2,                              &
                       klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,soot_aged)
  CASE (114)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(ocff_fresh(i1:i2,j1:j2,                             &
                        klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,ocff_fresh)
  CASE (115)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(ocff_aged(i1:i2,j1:j2,                              &
                       klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,ocff_aged)
  CASE (132)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(dms_sea_conc(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,dms_sea_conc)
  CASE (150)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(vertvel(i1:i2,j1:j2,0:ukcaD1codes(n)%len_dim3-1))
    CALL ukca_extract_d1_data(                                   &
     first,n,vertvel)
  CASE (211)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(conv_cloud_amount(i1:i2,j1:j2,                      &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,conv_cloud_amount)
  CASE (216)
    ALLOCATE(frac_types(ukcaD1codes(n)%len_dim1,                 &
                        ukcaD1codes(n)%len_dim2))
    CALL ukca_extract_d1_data(                                   &
     first,n,frac_types)
  CASE (217)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(laift_lp(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,laift_lp)
  CASE (218)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(canhtft_lp(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,canhtft_lp)
  CASE (229)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(canwctile_lp(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,canwctile_lp)
  CASE (233)
    ALLOCATE(tstar_tile(ukcaD1codes(n)%len_dim1,                 &
                        ukcaD1codes(n)%len_dim2))
    CALL ukca_extract_d1_data(                                   &
     first,n,tstar_tile)
  CASE (234)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(z0tile_lp(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,z0tile_lp)
  CASE (240)
    ALLOCATE(snow_tile(ukcaD1codes(n)%len_dim1,                  &
                       ukcaD1codes(n)%len_dim2))
    CALL ukca_extract_d1_data(                                   &
     first,n,snow_tile)
  CASE (252)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(co2_interactive(i1:i2,j1:j2,                        &
                        0:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,co2_interactive)
  CASE (253)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(rho_r2(i1:i2,j1:j2,                                 &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,rho_r2)
  CASE (254)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(qcl(i1:i2,j1:j2,                                    &
                        klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,qcl)
  CASE (255)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(exner_rho_levels(i1:i2,j1:j2,                       &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,exner_rho_levels)
  CASE (265)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(area_cloud_fraction(i1:i2,j1:j2,                    &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,area_cloud_fraction)
  CASE (266)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(cloud_frac(i1:i2,j1:j2,                             &
                        klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,cloud_frac)
  CASE (267)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(cloud_liq_frac(i1:i2,j1:j2,                         &
                            klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,cloud_liq_frac)
  CASE (351)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(biogenic(i1:i2,j1:j2,                               &
                      klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,biogenic)
  CASE (359)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(so4_accum(i1:i2,j1:j2,                              &
                        klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,so4_accum)
  CASE (360)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(so4_aitken(i1:i2,j1:j2,                             &
                        klev1:ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,so4_aitken)
  CASE (505)
    ALLOCATE(fland(ukcaD1codes(n)%len_dim1))
    CALL ukca_extract_d1_data(                                   &
     first,n,fland)
  CASE (510)
    ! Only fill radiation fields on radiation timesteps
    IF (l_rad_step_prog) THEN 
      CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
      IF ( .NOT. ALLOCATED(land_albedo))                         &
        ALLOCATE(land_albedo(i1:i2,j1:j2))
      CALL ukca_extract_d1_data(                                 &
      first,n,land_albedo)
    END IF

  CASE DEFAULT
    cmessage='Item not found in prognostic case statement'
    errcode = ukcaD1codes(n)%item
    CALL ereport('GETD1FLDS',errcode,cmessage)
  END SELECT

  ! Diagnostics (section 0): fill appropriate array
ELSE IF (UkcaD1codes(n)%section == 0 .AND.                       &
        .NOT. UkcaD1codes(n)%prognostic) THEN
  SELECT CASE(UkcaD1codes(n)%item)

  CASE (406)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(exner_theta_levels(i1:i2,j1:j2,                     &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,exner_theta_levels)
  CASE (407)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(p_rho_levels(i1:i2,j1:j2,                           &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,p_rho_levels)
  CASE (408)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(p_theta_levels(i1:i2,j1:j2,                         &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,p_theta_levels)
  CASE (409)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(pstar(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,pstar)
  CASE DEFAULT
    cmessage='N not found in diagnostic(0) case statement'
    errcode = n
    CALL ereport('GETD1FLDS',errcode,cmessage)
  END SELECT

  ! Diagnostics (section 1): fill appropriate array
ELSE IF (UkcaD1codes(n)%section == 1) THEN
  SELECT CASE(UkcaD1codes(n)%item)

  CASE (201)
    ! Only fill radiation fields on radiation timesteps
    IF (l_rad_step_prog) THEN 
      CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
      ALLOCATE(net_surf_SW(i1:i2,j1:j2))
      CALL ukca_extract_d1_data(                                   &
       first,n,net_surf_SW)
    END IF
  CASE (235)
    ! Only fill radiation fields on radiation timesteps
    IF (l_rad_step_prog) THEN 
      CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
      ALLOCATE(tot_surf_SW(i1:i2,j1:j2))
      CALL ukca_extract_d1_data(                                   &
       first,n,tot_surf_SW)
    END IF
  CASE (247)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(sea_salt_film(i1:i2,j1:j2,                            &
                           ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                     &
     first,n,sea_salt_film)
  CASE (248)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(sea_salt_jet(i1:i2,j1:j2,                             &
                          ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                     &
     first,n,sea_salt_jet)
  CASE DEFAULT
    cmessage='N not found in diagnostic(1) case statement'
    errcode = n
    CALL ereport('GETD1FLDS',errcode,cmessage)
  END SELECT

  ! Diagnostics (section 2): fill appropriate array
ELSE IF (UkcaD1codes(n)%section == 2) THEN
  SELECT CASE(UkcaD1codes(n)%item)

  CASE (284)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(sulphate_od(i1:i2,j1:j2,                            &
                       ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,sulphate_od)
  CASE DEFAULT
    cmessage='N not found in diagnostic(2) case statement'
    errcode = n
    CALL ereport('GETD1FLDS',errcode,cmessage)
  END SELECT

  !     Diagnostics (section 3): fill appropriate array

ELSE IF (UkcaD1codes(n)%section == 3) THEN
  SELECT CASE(UkcaD1codes(n)%item)
  CASE (60)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(rhokh_mix(i1:i2,j1:j2,                              &
                       ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,rhokh_mix)
  CASE (64)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(dtrdz_charney_grid(i1:i2,j1:j2,                     &
                       ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,dtrdz_charney_grid)
  CASE (65)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(kent(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,kent)
  CASE (66)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(we_lim(i1:i2,j1:j2,                                 &
                       ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,we_lim)
  CASE (67)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(t_frac(i1:i2,j1:j2,                                 &
                       ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,t_frac)
  CASE (68)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(zrzi(i1:i2,j1:j2,                                   &
                       ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,zrzi)
  CASE (69)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(kent_dsc(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,kent_dsc)
  CASE (70)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(we_lim_dsc(i1:i2,j1:j2,                             &
                       ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,we_lim_dsc)
  CASE (71)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(t_frac_dsc(i1:i2,j1:j2,                             &
                       ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,t_frac_dsc)
  CASE (72)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(zrzi_dsc(i1:i2,j1:j2,                               &
                       ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,zrzi_dsc)
  CASE (73)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(zhsc(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,zhsc)
  CASE (230)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(u_scalar_10m(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,u_scalar_10m)
  CASE (217)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(surf_hf(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,surf_hf)
  CASE (462)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(stcon(i1:i2,j1:j2,Ukcad1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,stcon)
  CASE (465)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(u_s(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,u_s)
  CASE (473)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(bl_tke(i1:i2,j1:j2,ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,bl_tke)
  CASE DEFAULT
    cmessage='N not found in diagnostic(3) case statement'
    errcode = n
    CALL ereport('GETD1FLDS',errcode,cmessage)
  END SELECT

  ! Diagnostics (section 4): fill appropriate array
ELSE IF (UkcaD1codes(n)%section == 4) THEN
  SELECT CASE(UkcaD1codes(n)%item)
  CASE (205)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(cloud_liq_water(i1:i2,j1:j2,                        &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,cloud_liq_water)
  CASE (222)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(ls_rain3d(i1:i2,j1:j2,                              &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,ls_rain3d)
  CASE (223)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(ls_snow3d(i1:i2,j1:j2,                              &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,ls_snow3d)
  CASE (227)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(ls_ppn_frac(i1:i2,j1:j2,                            &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,ls_ppn_frac)
  CASE (253)  
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)  
    ALLOCATE(ice_melt(i1:i2,j1:j2,                               &  
                        ukcaD1codes(n)%len_dim3))  
    CALL ukca_extract_d1_data(                                   &  
     first,n,ice_melt)  
  CASE (254)  
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)  
    ALLOCATE(snow_melt(i1:i2,j1:j2,                              &  
                        ukcaD1codes(n)%len_dim3))  
    CALL ukca_extract_d1_data(                                   &  
     first,n,snow_melt)  
  CASE (247)  
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)  
    ALLOCATE(rim_cry(i1:i2,j1:j2,                                &  
                      ukcaD1codes(n)%len_dim3))  
    CALL ukca_extract_d1_data(                                   &  
     first,n,rim_cry)  
  CASE (248)  
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)  
    ALLOCATE(rim_agg(i1:i2,j1:j2,                                &  
                      ukcaD1codes(n)%len_dim3))  
    CALL ukca_extract_d1_data(                                   &  
     first,n,rim_agg)  
  CASE (257)  
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)  
    ALLOCATE(autoconv(i1:i2,j1:j2,                               &  
                        ukcaD1codes(n)%len_dim3))  
    CALL ukca_extract_d1_data(                                   &  
     first,n,autoconv)  
  CASE (258)  
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)  
    ALLOCATE(accretion(i1:i2,j1:j2,                              &  
                        ukcaD1codes(n)%len_dim3))  
    CALL ukca_extract_d1_data(                                   &  
     first,n,accretion)
  CASE DEFAULT
    cmessage='N not found in diagnostic(4) case statement'
    errcode=1
    WRITE(umMessage,'(A72,I10)') cmessage,n
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    CALL ereport('GETD1FLDS',errcode,cmessage)
  END SELECT

  ! Diagnostics (section 5): fill appropriate array

ELSE IF (UkcaD1codes(n)%section == 5) THEN
  SELECT CASE(UkcaD1codes(n)%item)
  CASE (227)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(conv_rain3d(i1:i2,j1:j2,                            &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,conv_rain3d)
  CASE (228)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(conv_snow3d(i1:i2,j1:j2,                            &
                        ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,conv_snow3d)
  CASE (218)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(conv_cloud_base(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,conv_cloud_base)
  CASE (219)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(conv_cloud_top(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,conv_cloud_top)
  CASE DEFAULT
    cmessage='N not found in diagnostic(5) case statement'
    errcode=1

    CALL ereport('GETD1FLDS',errcode,cmessage)
  END SELECT

  ! Diagnostics (section 8): fill appropriate array

ELSE IF (UkcaD1codes(n)%section == 8) THEN
  SELECT CASE(UkcaD1codes(n)%item)
  CASE (242)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(ch4_wetl_emiss(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                   &
     first,n,ch4_wetl_emiss)
  CASE DEFAULT
    cmessage='N not found in diagnostic(8) case statement'
    errcode=1

    CALL ereport('GETD1FLDS',errcode,cmessage)
  END SELECT

  ! Diagnostics (section 15): fill appropriate array

ELSE IF (UkcaD1codes(n)%section == 15) THEN
  SELECT CASE(UkcaD1codes(n)%item)
  CASE (218)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(PV_on_theta_mlevs(i1:i2,j1:j2,                      &
             ukcaD1codes(n)%len_dim3))
    CALL ukca_extract_d1_data(                                   &
     first,n,PV_on_theta_mlevs)
  CASE DEFAULT
    cmessage='N not found in diagnostic(8) case statement'
    errcode=1

    CALL ereport('GETD1FLDS',errcode,cmessage)
  END SELECT

  ! Diagnostics (Section 30): fill appropriate array
ELSE IF (UkcaD1codes(n)%section == 30) THEN
  SELECT CASE(UkcaD1codes(n)%item)
  CASE (453)
    CALL ukca_set_array_bounds(n,i1,i2,j1,j2)
    ALLOCATE(tropopause_height(i1:i2,j1:j2))
    CALL ukca_extract_d1_data(                                    &
     first,n,tropopause_height)
  CASE DEFAULT
    cmessage='N not found in diagnostic(30) case statement'
    errcode = n
    CALL ereport('GETD1FLDS',errcode,cmessage)
  END SELECT

! If this ELSE is reached, then the program has failed to find 
! an array to allocate the data to in this IF block
ELSE
  cmessage='N not located in IF statement'
  errcode = n
  CALL ereport('GETD1FLDS',errcode,cmessage)
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE getd1flds

! ----------------------------------------------------------------------
! Put tracer and non-tracer prognostic fields back into D1 after UKCA
! ----------------------------------------------------------------------
SUBROUTINE putd1flds()
! ----------------------------------------------------------------------

USE ukca_set_array_bounds_mod, ONLY: ukca_set_array_bounds
USE ukca_option_mod, ONLY: l_ukca_chem
IMPLICIT NONE

INTEGER    :: n, i,i1,i2        ! loop variables
INTEGER    :: j,j1,j2           ! loop variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PUTD1FLDS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (L_ukca_chem) THEN

  ! Non transported prognostics

  DO i = 1, dim_ntp

    ! If ntp allocated, we need to copy it back to D1
    IF (ALLOCATED(all_ntp(i)%data_3d)) THEN 

      ! Look up location in ukcad1codes using section and item number 
      n = ukca_d1_myindex(all_ntp(i)%section , all_ntp(i)%item)

      ! If ENDGame - size of array in D1 = 0:model_levels i.e. 
      !  ntp_data(:,:,0:model_levels), but only 1:model_levels is
      !  filled as data_3d(:,:,1:model_levels), so copy data from 
      !  1st level to 0th level
      ntp_data(:,:,1:model_levels) =                              &
        all_ntp(i)%data_3d(:,:,1:model_levels)
      ntp_data(:,:,0) = ntp_data(:,:,1)

      ! put data into D1
      d1(ukcaD1codes(n)%address:ukcaD1codes(n)%address+           &
         ukcaD1codes(n)%length-1)=RESHAPE(ntp_data,(/ukcaD1codes(n)%length/))

   END IF
  END DO
  
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE putd1flds
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! To help with debugging, print values input from UM 
! ----------------------------------------------------------------------
SUBROUTINE ukca_pr_inputs(first)
! ----------------------------------------------------------------------
USE UM_ParCore,           ONLY: mype
USE umPrintMgr,           ONLY: umMessage, UmPrint, PrintStatus,  &
                                PrStatus_Diag, PrStatus_Oper
USE ukca_option_mod,      ONLY: l_ukca_dust,                      &
                                l_ukca_aerchem, l_ukca_nr_aqchem, &
                                l_ukca_offline_be, l_ukca_offline,&
                                l_ukca_classic_hetchem,           &
                                i_ukca_photol
                                
USE dyn_coriolis_mod,     ONLY: f3_at_u
USE rad_input_mod,        ONLY: l_use_biogenic,                   &
                                l_use_seasalt_direct,             & 
                                l_use_seasalt_indirect
USE mphys_inputs_mod,     ONLY: l_use_seasalt_autoconv
USE spec_sw_lw,           ONLY: sw_spectrum
USE run_aerosol_mod,      ONLY: l_sulpc_so2, l_soot, l_ocff,      &
                                l_use_seasalt_sulpc,              &
                                l_use_seasalt_sulpc, l_sulpc_so2
USE ukca_photo_scheme_mod,  ONLY: i_ukca_fastjx 
USE trignometric_mod,     ONLY: true_longitude,                   &
                                true_latitude,                    &
                                sin_theta_longitude,              &
                                sin_theta_latitude,               &
                                sin_v_latitude,                   &
                                cos_v_latitude,                   &
                                FV_cos_theta_latitude,            &
                                cos_theta_longitude,              &
                                tan_theta_latitude

IMPLICIT NONE

LOGICAL, INTENT(IN) :: first

INTEGER :: j, k ! loop counters
INTEGER :: klev1 

klev1   = tdims%k_start

! Debug these every TS only  IF (PrintStatus >= PrStatus_Diag)
IF (PrintStatus >= PrStatus_Diag) THEN
  
  ! print out info on NTP
  CALL print_all_ntp(all_ntp)

END IF
IF (first .AND. PrintStatus >= PrStatus_Oper ) THEN

  WRITE(umMessage,'(A)') ' ==================================='
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A)') '  MAX and MIN of UKCA INPUTS from D1'
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A)') ' ==================================='
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8, 2E12.4)') 'soil_moisture_layer1: ',mype,              &
    MAXVAL(soil_moisture_layer1),MINVAL(soil_moisture_layer1)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  !

  IF (L_ukca_dust) THEN
    WRITE(umMessage,'(A,I8, 2E12.4)') 'frac_types: ',mype,                     &
      MAXVAL(frac_types),MINVAL(frac_types)
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, 2E12.4)') 'tstar_tile: ',mype,                     &
      MAXVAL(tstar_tile),MINVAL(tstar_tile)
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, 2E12.4)') 'snow_tile: ',mype,                      &
      MAXVAL(snow_tile),MINVAL(snow_tile)
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  END IF

  WRITE(umMessage,'(A,I8, 2E12.4)') 'Rough_length: ',mype,                     &
    MAXVAL(Rough_length),MINVAL(Rough_length)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8, 2I8)') 'conv_cloud_base: ',mype,                     &
    MAXVAL(conv_cloud_base(:,:)),                                              &
    MINVAL(conv_cloud_base(:,:))
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8, 2I8)') 'conv_cloud_top: ',mype,                      &
    MAXVAL(conv_cloud_top(:,:)),                                               &
    MINVAL(conv_cloud_top(:,:))
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,L5)') 'land_sea_mask(1,1): ', land_sea_mask(1,1)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8, 2E12.4)') 'f3_at_u: ',mype,                          &
    MAXVAL(f3_at_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end)),      &
    MINVAL(f3_at_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end))
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8, 2E12.4)') 'FV_cos_theta_latitude: ',mype,            &
    MAXVAL(FV_cos_theta_latitude(:,:)),                                        &
    MINVAL(FV_cos_theta_latitude(:,:))
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')

  DO k=1,model_levels
    WRITE(umMessage,'(A,I6)') 'LEVEL: ',k
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'p_rho_levels:   ',mype, k,          &
      MAXVAL(p_rho_levels(1:row_length,1:rows,k)),                             &
      MINVAL(p_rho_levels(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'p_theta_levels: ',mype,k,           &
      MAXVAL(p_theta_levels(1:row_length,1:rows,k)),                           &
      MINVAL(p_theta_levels(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'rho_r2:         ',mype,k,           &
      MAXVAL(rho_r2(1:row_length,1:rows,k)),                                   &
      MINVAL(rho_r2(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    IF (ALLOCATED(co2_interactive)) THEN
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'co2_interactive: ',mype,k,        &
      MAXVAL(co2_interactive(1:row_length,1:rows,k)),                          &
      MINVAL(co2_interactive(1:row_length,1:rows,k))
      CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    END IF
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'ls_rain3d:      ',mype,k,           &
      MAXVAL(ls_rain3d(1:row_length,1:rows,k)),                                &
      MINVAL(ls_rain3d(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'ls_snow3d:      ',mype,k,           &
      MAXVAL(ls_snow3d(1:row_length,1:rows,k)),                                &
      MINVAL(ls_snow3d(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'conv_rain3d:      ',mype,k,         &
      MAXVAL(conv_rain3d(1:row_length,1:rows,k)),                              &
      MINVAL(conv_rain3d(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'conv_snow3d:      ',mype,k,         &
      MAXVAL(conv_snow3d(1:row_length,1:rows,k)),                              &
      MINVAL(conv_snow3d(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'PV_on_theta_mlevs:      ',mype,k,   &
      MAXVAL(PV_on_theta_mlevs(1:row_length,1:rows,k)),                        &
      MINVAL(PV_on_theta_mlevs(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')

    ! If heterogeneous chemistry on CLASSIC aerosols is used then
    ! print their corresponding MMRs or aerosol numbers.
    IF (l_ukca_classic_hetchem) THEN
        IF (l_soot) THEN
          WRITE(umMessage,'(A30,2I5,3E12.4)')                        &
                 'CLASSIC soot fresh MMR: ', mype, k,                &
                  MAXVAL(soot_fresh (1:row_length,1:rows,k)),        &
                  MINVAL(soot_fresh (1:row_length,1:rows,k)),        &
                  SUM   (soot_fresh (1:row_length,1:rows,k)) /       &
                  SIZE  (soot_fresh (1:row_length,1:rows,k))
          CALL umPrint(umMessage,src='ukca_main1-ukca_main1')

          WRITE(umMessage,'(A30,2I5,3E12.4)')                        &
                 'CLASSIC soot aged MMR: ', mype, k,                 &
                  MAXVAL(soot_aged (1:row_length,1:rows,k)),         &
                  MINVAL(soot_aged (1:row_length,1:rows,k)),         &
                  SUM   (soot_aged (1:row_length,1:rows,k)) /        &
                  SIZE  (soot_aged (1:row_length,1:rows,k))
          CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
        END IF

        IF (l_ocff) THEN
          WRITE(umMessage,'(A30,2I5,3E12.4)')                        &
                 'CLASSIC OCFF fresh MMR: ', mype, k,                &
                  MAXVAL(ocff_fresh (1:row_length,1:rows,k)),        &
                  MINVAL(ocff_fresh (1:row_length,1:rows,k)),        &
                  SUM   (ocff_fresh (1:row_length,1:rows,k)) /       &
                  SIZE  (ocff_fresh (1:row_length,1:rows,k))
          CALL umPrint(umMessage,src='ukca_main1-ukca_main1')

          WRITE(umMessage,'(A30,2I5,3E12.4)')                        &
                 'CLASSIC OCFF aged MMR: ', mype, k,                 &
                  MAXVAL(ocff_aged (1:row_length,1:rows,k)),         &
                  MINVAL(ocff_aged (1:row_length,1:rows,k)),         &
                  SUM   (ocff_aged (1:row_length,1:rows,k)) /        &
                  SIZE  (ocff_aged (1:row_length,1:rows,k))
          CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
        END IF

        IF (l_use_biogenic) THEN
          WRITE(umMessage,'(A30,2I5,3E12.4)')                        &
                 'CLASSIC biogenic aerosol MMR: ', mype, k,          &
                  MAXVAL(biogenic (1:row_length,1:rows,k)),          &
                  MINVAL(biogenic (1:row_length,1:rows,k)),          &
                  SUM   (biogenic (1:row_length,1:rows,k)) /         &
                  SIZE  (biogenic (1:row_length,1:rows,k))
          CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
        END IF

        IF (l_use_seasalt_autoconv .OR. l_use_seasalt_sulpc .OR.     &
            l_use_seasalt_indirect .OR. l_use_seasalt_direct) THEN
          WRITE(umMessage,'(A30,2I5,3E12.4)')                        &
                 'CLASSIC sea salt film no.: ', mype, k,             &
                  MAXVAL(sea_salt_film (1:row_length,1:rows,k)),     &
                  MINVAL(sea_salt_film (1:row_length,1:rows,k)),     &
                  SUM   (sea_salt_film (1:row_length,1:rows,k)) /    &
                  SIZE  (sea_salt_film (1:row_length,1:rows,k))
          CALL umPrint(umMessage,src='ukca_main1-ukca_main1')

          WRITE(umMessage,'(A30,2I5,3E12.4)')                        &
                 'CLASSIC sea salt jet no.: ', mype, k,              &
                  MAXVAL(sea_salt_jet (1:row_length,1:rows,k)),      &
                  MINVAL(sea_salt_jet (1:row_length,1:rows,k)),      &
                  SUM   (sea_salt_jet (1:row_length,1:rows,k)) /     &
                  SIZE  (sea_salt_jet (1:row_length,1:rows,k))
          CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
        END IF
    END IF

  END DO   ! k

  ! wet variables
  DO k=1,model_levels
    WRITE(umMessage,'(A,I6)') 'WET LEVEL: ',k
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'q:              ',mype,k,           &
      MAXVAL(q(1:row_length,1:rows,k)),                                        &
      MINVAL(q(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'qcl:            ',mype,k,           &
      MAXVAL(qcl(1:row_length,1:rows,k)),                                      &
      MINVAL(qcl(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'qcf:            ',mype,k,           &
      MAXVAL(qcf(1:row_length,1:rows,k)),                                      &
      MINVAL(qcf(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A12,2I5,4E12.4)') 'autoconv: ',mype,k,                   &
      MAXVAL(autoconv(1:row_length,1:rows,k)),                                 &
      MINVAL(autoconv(1:row_length,1:rows,k)),                                 &
      SUM(autoconv(1:row_length,1:rows,k))/                                    &
      SIZE(autoconv(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A12,2I5,4E12.4)') 'accretion: ',mype,k,                  &
      MAXVAL(accretion(1:row_length,1:rows,k)),                                &
      MINVAL(accretion(1:row_length,1:rows,k)),                                &
      SUM(accretion(1:row_length,1:rows,k))/                                   &
      SIZE(accretion(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A12,2I5,4E12.4)') 'rim_agg: ',mype,k,                    &
      MAXVAL(rim_agg(1:row_length,1:rows,k)),                                  &
      MINVAL(rim_agg(1:row_length,1:rows,k)),                                  &
      SUM(rim_agg(1:row_length,1:rows,k))/                                     &
      SIZE(rim_agg(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A12,2I5,4E12.4)') 'rim_cry: ',mype,k,                    &
      MAXVAL(rim_cry(1:row_length,1:rows,k)),                                  &
      MINVAL(rim_cry(1:row_length,1:rows,k)),                                  &
      SUM(rim_cry(1:row_length,1:rows,k))/                                     &
      SIZE(rim_cry(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A12,2I5,4E12.4)') 'ls_rain3d: ',mype,k,                  &
      MAXVAL(ls_rain3d(1:row_length,1:rows,k)),                                &
      MINVAL(ls_rain3d(1:row_length,1:rows,k)),                                &
      SUM(ls_rain3d(1:row_length,1:rows,k))/                                   &
      SIZE(ls_rain3d(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A12,2I5,4E12.4)') 'ls_snow : ',mype,k,                   &
      MAXVAL(ls_snow3d(1:row_length,1:rows,k)),                                &
      MINVAL(ls_snow3d(1:row_length,1:rows,k)),                                &
      SUM(ls_snow3d(1:row_length,1:rows,k))/                                   &
      SIZE(ls_snow3d(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  END DO

  ! Exner - additional level at top
  DO k=1,model_levels+1
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'EXNER_rho_levels: ',mype,k,         &
      MAXVAL(exner_rho_levels(1:row_length,1:rows,k)),                         &
      MINVAL(exner_rho_levels(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  END DO
  
  ! boundary layer variables  
  DO k=1,bl_levels
    WRITE(umMessage,'(A,I6)') 'BL LEVEL: ',k
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'rhokh_mix:     ',mype,k,            &
      MAXVAL(rhokh_mix(1:row_length,1:rows,k)),                                &
      MINVAL(rhokh_mix(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'dtrdz_charney_grid:     ',mype,k,   &
      MAXVAL(dtrdz_charney_grid(1:row_length,1:rows,k)),                       &
      MINVAL(dtrdz_charney_grid(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  END DO
  
  ! Sulfate aerosol inputs to photolysis schemes
  IF ( l_sulpc_so2 .AND. i_ukca_photol == i_ukca_fastjx ) THEN
    DO k=1,model_levels
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'so4_aitken: ', mype,k,            &
              MAXVAL(so4_aitken(1:row_length,1:rows,k)),                       &
              MINVAL(so4_aitken(1:row_length,1:rows,k))
      CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'so4_accum: ', mype,k,             &
              MAXVAL(so4_accum(1:row_length,1:rows,k)),                        &
              MINVAL(so4_accum(1:row_length,1:rows,k))
      CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    END DO
    IF ((i_ukca_photol == i_ukca_fastjx) .AND.                                 &
      ALLOCATED(sulphate_od)) THEN
      DO k=1,sw_spectrum(1)%basic%n_band
        WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'sulphate_od: ', mype,k,         &
          MAXVAL(sulphate_od(1:row_length,1:rows,k)),                          &
          MINVAL(sulphate_od(1:row_length,1:rows,k))
        CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
      END DO
    END IF
  END IF        ! sulpc_so2 and Fast-J/X

  WRITE(umMessage,'(A,I8,2I8)') 'kent:     ', mype,                            &
    MAXVAL(kent), MINVAL(kent)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8,2I8)') 'kent_dsc: ',mype,                             &
    MAXVAL(kent_dsc), MINVAL(kent_dsc)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8,2E12.4)') 'zhsc:     ', mype,                         &
    MAXVAL(zhsc), MINVAL(zhsc)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8,2E12.4)') 'BL_depth: ',mype,                          &
    MAXVAL(zbl), MINVAL(zbl)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  DO k=1,3
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'we_lim: ', mype,k,                  &
      MAXVAL(we_lim(1:row_length,1:rows,k)),                                   &
      MINVAL(we_lim(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 't_frac: ', mype,k,                  &
      MAXVAL(t_frac(1:row_length,1:rows,k)),                                   &
      MINVAL(t_frac(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'zrzi: ', mype,k,                    &
      MAXVAL(zrzi(1:row_length,1:rows,k)),                                     &
      MINVAL(zrzi(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'we_lim_dsc: ', mype,k,              &
      MAXVAL(we_lim_dsc(1:row_length,1:rows,k)),                               &
      MINVAL(we_lim_dsc(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 't_frac_dsc: ', mype,k,              &
      MAXVAL(t_frac_dsc(1:row_length,1:rows,k)),                               &
      MINVAL(t_frac_dsc(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'zrzi_dsc: ', mype,k,                &
      MAXVAL(zrzi_dsc(1:row_length,1:rows,k)),                                 &
      MINVAL(zrzi_dsc(1:row_length,1:rows,k))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  END DO
  WRITE(umMessage,'(A,I8,2E12.4)') 'pstar: ', mype,MAXVAL(pstar), MINVAL(pstar)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8,2E12.4)') 'Tstar: ', mype,MAXVAL(Tstar), MINVAL(Tstar)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8,2E12.4)') 'fland: ', mype,MAXVAL(fland), MINVAL(fland)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A,I8,2E12.4)') 'u_s  : ', mype,MAXVAL(u_s  ), MINVAL(u_s  )
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  WRITE(umMessage,'(A16,I5,2E12.4)') 'seaice_frac: ', mype,                    &
    MAXVAL(seaice_frac), MINVAL(seaice_frac)
  CALL umPrint(umMessage,src='ukca_main1-ukca_main1')

END IF 
! ----------------------------------------------------------------------
! End of debugging write statements
! ----------------------------------------------------------------------
END SUBROUTINE ukca_pr_inputs


!-------------------------------------------------------------------------------
SUBROUTINE ukca_um_dealloc()

USE ukca_all_tracers_copy_mod,  ONLY: ukca_q_increment
IMPLICIT NONE

! POSSIBLE THINGS TO DEALLOCATE
IF (ALLOCATED(sea_salt_jet))  DEALLOCATE(sea_salt_jet)
IF (ALLOCATED(sea_salt_film)) DEALLOCATE(sea_salt_film)
IF (ALLOCATED(biogenic))      DEALLOCATE(biogenic)
IF (ALLOCATED(ocff_aged))     DEALLOCATE(ocff_aged)
IF (ALLOCATED(ocff_fresh))    DEALLOCATE(ocff_fresh)
IF (ALLOCATED(soot_aged))     DEALLOCATE(soot_aged)
IF (ALLOCATED(soot_fresh))    DEALLOCATE(soot_fresh)
IF (ALLOCATED(net_surf_sw)) DEALLOCATE(net_surf_sw)
IF (ALLOCATED(tot_surf_sw)) DEALLOCATE(tot_surf_sw)
IF (ALLOCATED(trmol_post_atmstep)) DEALLOCATE(trmol_post_atmstep)  ! moles

! from PUTD1FLDS
IF (ALLOCATED(emission1)) DEALLOCATE(emission1)
IF (ALLOCATED(ntp_data)) DEALLOCATE(ntp_data)
IF (ALLOCATED(all_emissions)) DEALLOCATE(all_emissions)

! ACTUAL THINGS TO DEALLOCATE
IF (ALLOCATED(totnodens)) DEALLOCATE(totnodens)
IF (ALLOCATED(conv_cloud_base)) DEALLOCATE(conv_cloud_base)
IF (ALLOCATED(conv_cloud_top)) DEALLOCATE(conv_cloud_top)
IF (ALLOCATED(land_sea_mask)) DEALLOCATE(land_sea_mask)
IF (ALLOCATED(theta)) DEALLOCATE(theta)
IF (ALLOCATED(p_rho_levels)) DEALLOCATE(p_rho_levels)
IF (ALLOCATED(p_theta_levels)) DEALLOCATE(p_theta_levels)
IF (ALLOCATED(exner_theta_levels)) DEALLOCATE(exner_theta_levels)
IF (ALLOCATED(exner_rho_levels)) DEALLOCATE(exner_rho_levels)
IF (ALLOCATED(q)) DEALLOCATE(q)
IF (ALLOCATED(rho_r2)) DEALLOCATE(rho_r2)
IF (ALLOCATED(co2_interactive)) DEALLOCATE(co2_interactive)
IF (ALLOCATED(qcl)) DEALLOCATE(qcl)
IF (ALLOCATED(qcf)) DEALLOCATE(qcf)
IF (ALLOCATED(ls_ppn_frac)) DEALLOCATE(ls_ppn_frac)
IF (ALLOCATED(cloud_liq_water)) DEALLOCATE(cloud_liq_water)
IF (ALLOCATED(ls_rain3d)) DEALLOCATE(ls_rain3d)
IF (ALLOCATED(ls_snow3d)) DEALLOCATE(ls_snow3d)
IF (ALLOCATED(autoconv)) DEALLOCATE(autoconv)  
IF (ALLOCATED(accretion)) DEALLOCATE(accretion)
IF (ALLOCATED(snow_melt)) DEALLOCATE(snow_melt) 
IF (ALLOCATED(ice_melt)) DEALLOCATE(ice_melt)
IF (ALLOCATED(rim_cry)) DEALLOCATE(rim_cry) 
IF (ALLOCATED(rim_agg)) DEALLOCATE(rim_agg)
IF (ALLOCATED(conv_rain3d)) DEALLOCATE(conv_rain3d)
IF (ALLOCATED(conv_snow3d)) DEALLOCATE(conv_snow3d)
IF (ALLOCATED(conv_cloud_amount)) DEALLOCATE(conv_cloud_amount)
IF (ALLOCATED(area_cloud_fraction)) DEALLOCATE(area_cloud_fraction)
IF (ALLOCATED(so4_aitken)) DEALLOCATE(so4_aitken)
IF (ALLOCATED(so4_accum)) DEALLOCATE(so4_accum)
IF (ALLOCATED(sulphate_od)) DEALLOCATE(sulphate_od)
IF (ALLOCATED(dms_sea_conc)) DEALLOCATE(dms_sea_conc)
IF (ALLOCATED(chloro_sea)) DEALLOCATE(chloro_sea)
IF (ALLOCATED(rel_humid_frac)) DEALLOCATE(rel_humid_frac)
IF (ALLOCATED(qsvp)) DEALLOCATE(qsvp)
IF (ALLOCATED(U_scalar_10m)) DEALLOCATE(U_scalar_10m)
IF (ALLOCATED(pstar)) DEALLOCATE(pstar)
IF (ALLOCATED(tstar)) DEALLOCATE(tstar)
IF (ALLOCATED(dust_flux)) DEALLOCATE(dust_flux)
IF (ALLOCATED(soil_moisture_layer1)) DEALLOCATE(soil_moisture_layer1)
IF (ALLOCATED(Tstar_tile)) DEALLOCATE(Tstar_tile)
IF (ALLOCATED(Snow_tile)) DEALLOCATE(Snow_tile)
IF (ALLOCATED(Frac_types)) DEALLOCATE(Frac_types)
IF (ALLOCATED(Rough_length)) DEALLOCATE(Rough_length)
IF (ALLOCATED(fland)) DEALLOCATE(fland)
IF (ALLOCATED(zbl)) DEALLOCATE(zbl)
IF (ALLOCATED(surf_hf)) DEALLOCATE(surf_hf)
IF (ALLOCATED(seaice_frac)) DEALLOCATE(seaice_frac)
IF (ALLOCATED(conv_cloud_lwp)) DEALLOCATE(conv_cloud_lwp)
IF (ALLOCATED(u_s)) DEALLOCATE(u_s)
IF (ALLOCATED(stcon)) DEALLOCATE(stcon)
IF (ALLOCATED(laift_lp)) DEALLOCATE(laift_lp)
IF (ALLOCATED(canhtft_lp)) DEALLOCATE(canhtft_lp)
IF (ALLOCATED(z0tile_lp)) DEALLOCATE(z0tile_lp)
IF (ALLOCATED(canwctile_lp)) DEALLOCATE(canwctile_lp)
IF (ALLOCATED(ch4_wetl_emiss)) DEALLOCATE(ch4_wetl_emiss)
IF (ALLOCATED(tropopause_height)) DEALLOCATE(tropopause_height)
IF (ALLOCATED(strat_fluxdiags)) DEALLOCATE(strat_fluxdiags)
IF (ALLOCATED(mode_diags)) DEALLOCATE(mode_diags)
IF (ALLOCATED(cloud_frac)) DEALLOCATE(cloud_frac)
IF (ALLOCATED(dtrdz_charney_grid)) DEALLOCATE(dtrdz_charney_grid)
IF (ALLOCATED(kent)) DEALLOCATE(kent)
IF (ALLOCATED(kent_dsc)) DEALLOCATE(kent_dsc)
IF (ALLOCATED(PV_on_theta_mlevs)) DEALLOCATE(PV_on_theta_mlevs)
IF (ALLOCATED(rhokh_mix)) DEALLOCATE(rhokh_mix)
IF (ALLOCATED(bl_tke)) DEALLOCATE(bl_tke)
IF (ALLOCATED(vertvel)) DEALLOCATE(vertvel)
IF (ALLOCATED(t_frac)) DEALLOCATE(t_frac)
IF (ALLOCATED(t_frac_dsc)) DEALLOCATE(t_frac_dsc)
IF (ALLOCATED(um_ozone)) DEALLOCATE(um_ozone)
IF (ALLOCATED(we_lim)) DEALLOCATE(we_lim)
IF (ALLOCATED(we_lim_dsc)) DEALLOCATE(we_lim_dsc)
IF (ALLOCATED(zhsc)) DEALLOCATE(zhsc)
IF (ALLOCATED(zrzi)) DEALLOCATE(zrzi)
IF (ALLOCATED(zrzi_dsc)) DEALLOCATE(zrzi_dsc)
IF (ALLOCATED(cloud_liq_frac)) DEALLOCATE(cloud_liq_frac)
IF (ALLOCATED(ukca_q_increment)) DEALLOCATE(ukca_q_increment)
CALL ntp_dealloc(all_ntp)

END SUBROUTINE ukca_um_dealloc

!-------------------------------------------------------------------------------
END MODULE ukca_um_interf_mod
