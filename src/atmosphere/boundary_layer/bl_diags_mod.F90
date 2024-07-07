! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  ---------------------------------------------------------------------
!  Sructure containing new boundary layer diagnostics.
!  This permits easier addition of new boundary layer
!  diagnostics without additional passing of arguments
!  though the boundary layer tree.
!  It also does not require the addition
!  of extra subroutine arguments when adding a new diagnostic.

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer

!- ----------------------------------------------------------------------

MODULE bl_diags_mod

IMPLICIT NONE
SAVE

TYPE strnewbldiag

  ! Need to create a flag and a pointer

  LOGICAL :: l_ftl
  LOGICAL :: l_fqw
  LOGICAL :: l_taux
  LOGICAL :: l_tauy
  LOGICAL :: l_t_incr
  LOGICAL :: l_q_incr
  LOGICAL :: l_qcl_incr
  LOGICAL :: l_qcf_incr
  LOGICAL :: l_u_incr
  LOGICAL :: l_v_incr
  LOGICAL :: l_w_incr
  LOGICAL :: l_dtfric
  LOGICAL :: l_tl_incr
  LOGICAL :: l_qtl_incr
  LOGICAL :: l_cf_incr
  LOGICAL :: l_cfl_incr
  LOGICAL :: l_cff_incr
  LOGICAL :: l_smltop
  LOGICAL :: l_dsctop
  LOGICAL :: l_zhlocal
  LOGICAL :: l_zhpar
  LOGICAL :: l_dscbase
  LOGICAL :: l_cldbase
  LOGICAL :: l_weparm
  LOGICAL :: l_weparm_dsc
  LOGICAL :: l_dzh
  LOGICAL :: l_oblen
  LOGICAL :: l_ustar
  LOGICAL :: l_wbsurf
  LOGICAL :: l_gradrich
  LOGICAL :: l_wstar
  LOGICAL :: l_dbdz
  LOGICAL :: l_dvdzm
  LOGICAL :: l_rhokm
  LOGICAL :: l_rhokh
  LOGICAL :: l_tke
  LOGICAL :: l_ostressx
  LOGICAL :: l_ostressy
  LOGICAL :: l_elm3d
  LOGICAL :: l_elh3d
  LOGICAL :: l_rhokmloc
  LOGICAL :: l_rhokhloc
  LOGICAL :: l_rhokmsurf
  LOGICAL :: l_rhokhsurf
  LOGICAL :: l_rhokmsc
  LOGICAL :: l_rhokhsc
  LOGICAL :: l_weight1d
  LOGICAL :: l_fm
  LOGICAL :: l_fh
  LOGICAL :: l_rhogamu
  LOGICAL :: l_rhogamv
  LOGICAL :: l_rhogamt
  LOGICAL :: l_rhogamq
  LOGICAL :: l_elm
  LOGICAL :: l_tke_shr_prod
  LOGICAL :: l_tke_boy_prod
  LOGICAL :: l_tke_dissp
  LOGICAL :: l_sm
  LOGICAL :: l_sh
  LOGICAL :: l_wb_ng
  LOGICAL :: l_cf_trb
  LOGICAL :: l_ql_trb
  LOGICAL :: l_sgm_trb

  REAL, ALLOCATABLE :: t_incr(:,:,:)
  !                    181      temperature increment
  REAL, ALLOCATABLE :: q_incr(:,:,:)
  !                    182      vapour increment
  REAL, ALLOCATABLE :: qcl_incr(:,:,:)
  !                    183      liquid water increment
  REAL, ALLOCATABLE :: qcf_incr(:,:,:)
  !                    184      ice water increment
  REAL, ALLOCATABLE :: u_incr(:,:,:)
  !                    185      u wind increment
  REAL, ALLOCATABLE :: v_incr(:,:,:)
  !                    186      v wind increment
  REAL, ALLOCATABLE :: w_incr(:,:,:)
  !                    187      w wind increment
  REAL, ALLOCATABLE :: dTfric(:, :, :)
  !                    188      Heating increment from turbulence dissipation
  REAL, ALLOCATABLE :: cf_incr(:,:,:)
  !                    192      bulk cloud fraction increment
  REAL, ALLOCATABLE :: cfl_incr(:,:,:)
  !                    193      liquid cloud fraction increment
  REAL, ALLOCATABLE :: cff_incr(:,:,:)
  !                    194      frozen cloud fraction increment
  REAL, ALLOCATABLE :: smltop(:, :)
  !                    356      Top of surface mixed layer
  REAL, ALLOCATABLE :: dsctop(:, :)
  !                    357      Top of decoupled stratocu layer
  REAL, ALLOCATABLE :: zhlocal(:, :)
  !                    358      BL depth diagnosed from Ri>RiCrit
  REAL, ALLOCATABLE :: zhpar(:, :)
  !                    359      Height of diagnosis parcel top
  REAL, ALLOCATABLE :: dscbase(:, :)
  !                    360      Height of decoupled layer base
  REAL, ALLOCATABLE :: cldbase(:, :)
  !                    361      Height of stratocumulus cloud base
  REAL, ALLOCATABLE :: weparm(:, :)
  !                    362      Entrainment rate for SML
  REAL, ALLOCATABLE :: weparm_dsc(:, :)
  !                    363      Entrainment rate for DSC
  REAL, ALLOCATABLE :: dzh(:, :)
  !                    364      Inversion thickness
  REAL, ALLOCATABLE :: oblen(:, :)
  !                    464      Surface Obukhov length
  REAL, ALLOCATABLE :: ustar(:, :)
  !                    465      Friction velocity
  REAL, ALLOCATABLE :: wstar(:, :)
  !                    466      Convective velocity scale
  REAL, ALLOCATABLE :: wbsurf(:, :)
  !                    467      Surface buoyancy flux
  REAL, ALLOCATABLE :: gradrich(:, :, :)
  !                    468      Gradient Richardson number
  REAL, ALLOCATABLE :: dbdz(:, :, :)
  !                    469      Vertical buoyancy gradient
  REAL, ALLOCATABLE :: dvdzm(:, :, :)
  !                    470      Modulus of wind shear
  REAL, ALLOCATABLE :: rhokm(:, :, :)
  !                    471      BL Momentum diffusivity
  REAL, ALLOCATABLE :: rhokh(:, :, :)
  !                    472      BL Thermal diffusivity
  REAL, ALLOCATABLE :: tke(:, :, :)
  !                    473      Turbulent kinetic energy
  REAL, ALLOCATABLE :: ostressx(:, :, :)
  !                    474      Orographic stress (x-component)
  REAL, ALLOCATABLE :: ostressy(:, :, :)
  !                    475      Orographic stress (y-component)
  REAL, ALLOCATABLE :: elm3d(:, :, :)
  !                    501      Mixing length for momentum
  REAL, ALLOCATABLE :: elh3d(:, :, :)
  !                    502      Mixing length for heat and moisture
  REAL, ALLOCATABLE :: rhokmloc(:, :, :)
  !                    503      Km diffusion coeff from local scheme
  REAL, ALLOCATABLE :: rhokhloc(:, :, :)
  !                    504      Kh diffusion coeff from local scheme
  REAL, ALLOCATABLE :: rhokmsurf(:, :, :)
  !                    505      Km diffusion coeff for surface-driven turb
  REAL, ALLOCATABLE :: rhokhsurf(:, :, :)
  !                    506      Kh diffusion coeff for surface-driven turb
  REAL, ALLOCATABLE :: rhokmsc(:, :, :)
  !                    507      Km diffusion coeff for  cloud-top-driven turb
  REAL, ALLOCATABLE :: rhokhsc(:, :, :)
  !                    508      Kh diffusion coeff for  cloud-top-driven turb
  REAL, ALLOCATABLE :: fh(:, :, :)
  !                    511      stability function for scalars
  REAL, ALLOCATABLE :: fm(:, :, :)
  !                    512      stability function for momentum
  REAL, ALLOCATABLE :: weight1d(:, :, :)
  !                    513      weighting applied to 1D BL scheme in Smag
  !                             blending
  REAL, ALLOCATABLE :: rhogamu(:, :, :)
  !                    130      counter gradient term of taux
  REAL, ALLOCATABLE :: rhogamv(:, :, :)
  !                    131      counter gradient term of tauy
  REAL, ALLOCATABLE :: rhogamt(:, :, :)
  !                    132      counter gradient term of ftl
  REAL, ALLOCATABLE :: rhogamq(:, :, :)
  !                    133      counter gradient term of fqw
  REAL, ALLOCATABLE :: elm(:, :, :)
  !                    134      mixing length
  REAL, ALLOCATABLE :: tke_shr_prod(:, :, :)
  !                    135      production rate of TKE by shear
  REAL, ALLOCATABLE :: tke_boy_prod(:, :, :)
  !                    136      production rate of TKE by buoyancy
  REAL, ALLOCATABLE :: tke_dissp(:, :, :)
  !                    137      dissipation rate of TKE
  REAL, ALLOCATABLE :: sm(:, :, :)
  !                    138      non-dimensional diffusion coef. for u, v
  REAL, ALLOCATABLE :: sh(:, :, :)
  !                    139      non-dimensional diffusion coef. for t, q
  REAL, ALLOCATABLE :: wb_ng(:, :, :)
  !                    140      non-gradient buoyancy flux
  REAL, ALLOCATABLE :: cf_trb(:, :, :)
  !                    141      cloud fraction used in the TKE schemes
  REAL, ALLOCATABLE :: ql_trb(:, :, :)
  !                    142      condensed water used in the TKE schemes
  REAL, ALLOCATABLE :: sgm_trb(:, :, :)
  !                    143      standard deviation of the distribution
  !                             function in the TKE schemes

END TYPE strnewbldiag

TYPE (Strnewbldiag) :: BL_diag
! ----------------------------------------------------------------------
END MODULE bl_diags_mod
