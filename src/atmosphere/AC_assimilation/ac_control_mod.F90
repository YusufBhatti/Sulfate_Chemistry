! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation

! Description:
!   Declares variables which control the execution of the ac assimilation.

!   These variables previously declared in include files 
!   comacp, comag, comacdg, mppac

MODULE ac_control_mod

USE atmos_max_sizes, ONLY: model_levels_max, rows_max
USE comobs_mod, ONLY: nobtypmx

IMPLICIT NONE

! comacp

INTEGER,PARAMETER :: modeacp = 36
INTEGER,PARAMETER :: nanaltyp = 30
INTEGER,PARAMETER :: nradars = 15

INTEGER :: nact
INTEGER :: nprog
INTEGER :: ac_obs_types(nobtypmx)
INTEGER :: lact(nobtypmx)
INTEGER :: group_index(nobtypmx)
INTEGER :: type_index(nobtypmx)
INTEGER :: group_first(nobtypmx)
INTEGER :: group_last(nobtypmx)
INTEGER :: obs_unitno
INTEGER :: obs_format
INTEGER :: no_obs_files
INTEGER :: diag_rdobs
INTEGER :: iunitno
INTEGER :: mgeowt
INTEGER :: n_groups
INTEGER :: group_no(nobtypmx)
INTEGER :: mhcorfn
INTEGER :: macdiag(modeacp)
INTEGER :: mwtfn
INTEGER :: mdatadfn
INTEGER :: npass_rf
INTEGER :: nslabs_scfact(model_levels_max)
INTEGER :: no_scfact(nobtypmx)
INTEGER :: iomitobs(nanaltyp)
INTEGER :: master_ac_types(nobtypmx)
INTEGER :: def_ac_order(nobtypmx)
INTEGER :: def_no_iterations(nobtypmx)
INTEGER :: def_interval_iter(nobtypmx)
INTEGER :: def_no_anal_levs(nobtypmx)
INTEGER :: def_no_wt_levs(nobtypmx)
INTEGER :: def_mode_hanal(nobtypmx)
INTEGER :: lenact(nobtypmx)
INTEGER :: def_obthin(nobtypmx)
INTEGER :: mvint205
INTEGER :: mrampfn
INTEGER :: mglossfn
INTEGER :: lhn_range
INTEGER :: npass_rf_lhn
INTEGER :: wb_lonoffset
INTEGER :: wb_lonpts
INTEGER :: wb_latoffset
INTEGER :: wb_latpts

REAL :: obtime_nom
REAL :: vert_filt
REAL :: geowt_h(rows_max -1)
REAL :: troplat
REAL :: geowt_v(model_levels_max)
REAL :: vert_cor_scale(model_levels_max, 4)
REAL :: vert_cutoff_sl
REAL :: vert_cutoff_bw
REAL :: vert_cutoff_bh
REAL :: non_div_cor
REAL :: non_div_cor_10m
REAL :: speed_limit305
REAL :: tropint
REAL :: timef_start
REAL :: timef_obtime
REAL :: timef_end
REAL :: cscfact_h(rows_max)
REAL :: cscfact_v(model_levels_max)
REAL :: def_timeb(nobtypmx)
REAL :: def_timea(nobtypmx)
REAL :: def_tgetobb(nobtypmx)
REAL :: def_tgetoba(nobtypmx)
REAL :: def_cscale_start(nobtypmx)
REAL :: def_cscale_obtime(nobtypmx)
REAL :: def_cscale_end(nobtypmx)
REAL :: def_radinf(nobtypmx)
REAL :: wb_lat_cc(rows_max)
REAL :: wb_vert_v(model_levels_max)
REAL :: wb_land_factor
REAL :: radar_lat(nradars)
REAL :: radar_lon(nradars)
REAL :: radar_range_max
REAL :: epsilon_lhn
REAL :: relax_cf_lhn
REAL :: f1_506 , f2_506 , f3_506
REAL :: alpha_lhn
REAL :: lhn_limit
REAL :: fi_scale_lhn

REAL :: def_nudge_nh(nobtypmx)
REAL :: def_nudge_tr(nobtypmx)
REAL :: def_nudge_sh(nobtypmx)
REAL :: def_nudge_lam(nobtypmx)

REAL :: def_fi_var_factor(nobtypmx)
REAL :: fi_scale
REAL :: fi_scale_factor(model_levels_max)
REAL :: df_scale
REAL :: df_scale_lev(model_levels_max)
REAL :: df_coeff(model_levels_max)
REAL :: thresh_dl
REAL :: thresh_lm
REAL :: thresh_mh
REAL :: thresh_rmsf
REAL :: radar_range
REAL :: northlat, southlat, westlon, eastlon
REAL :: vert_cor_aero

LOGICAL :: lgeo
LOGICAL :: lhydr
LOGICAL :: lhydrol
LOGICAL :: lsyn
LOGICAL :: ltimer_ac
LOGICAL :: lac_uars
LOGICAL :: lac_mes
LOGICAL :: lwbal_sf,     lwbal_ua
LOGICAL :: wb_theta_ua, wb_land_scale, wb_theta_sf
LOGICAL :: lradar (nradars)
LOGICAL :: l_latlon_prver
LOGICAL :: l_mops_equals_rh
LOGICAL :: lcheck_grid
LOGICAL :: l_506_oberr
LOGICAL :: l_lhn , l_lhn_scale
LOGICAL :: l_lhn_1A
LOGICAL :: l_lhn_search , lhn_diag
LOGICAL :: l_verif_range
LOGICAL :: l_lhn_limit
LOGICAL :: l_lhn_fact
LOGICAL :: l_lhn_filt
LOGICAL :: l_obs_check
LOGICAL :: remove_neg_lh
LOGICAL :: use_conv_in_mops

! comag

INTEGER :: def_agres_rows(nobtypmx)
INTEGER :: def_agres_pts(nobtypmx)
INTEGER :: nrowsag
INTEGER :: nptsagmx
INTEGER :: nptsagmn
INTEGER :: nptsag(rows_max)
INTEGER :: npts0ag(rows_max+1)
INTEGER :: min_agpts

REAL :: stagrow1
REAL :: stagpt1
REAL :: row1mg
REAL :: row1mgth
REAL :: row1mguv
REAL :: row1ag
REAL :: dlatag
REAL :: dlatmg
REAL :: pt1mgth
REAL :: pt1mguv
REAL :: pt1ag
REAL :: dlongmg
REAL :: aglatdec
REAL :: agrowlen
REAL :: dlongag(rows_max)
REAL :: cosrowag(rows_max)

REAL :: dlat,dlong,xlatn,xlongw
REAL :: elfplat,elfplon

LOGICAL :: lagnp
LOGICAL :: lagsp

! comacdg

INTEGER,PARAMETER :: nldacp=5
INTEGER,PARAMETER :: ndacop=200
INTEGER,PARAMETER :: ndacp=100
INTEGER,PARAMETER :: ndacvp=11

LOGICAL :: ldiagac
LOGICAL :: lldac(nldacp)
LOGICAL :: lldag0
LOGICAL :: lldag(nldacp)
LOGICAL :: llband
LOGICAL :: ltitle
LOGICAL :: lrms
LOGICAL :: lnormf
LOGICAL :: ltemp
LOGICAL :: lverif

INTEGER :: mdiagtps(nobtypmx)
INTEGER :: mdiag
INTEGER :: ndgprt
INTEGER :: ndaco
INTEGER :: mdaco(ndacop)
INTEGER :: ndac
INTEGER :: ndacprt
INTEGER :: mdac(ndacp)
INTEGER :: ndacv
INTEGER :: modaco

REAL :: dlatn
REAL :: dlats
REAL :: dlongw
REAL :: dlonge
REAL :: daglat
REAL :: daglon
REAL :: dacv(ndacp,ndacvp)

CHARACTER(LEN=10) :: cacv(ndacvp)
CHARACTER(LEN=12) :: cact1
CHARACTER(LEN=32) :: cact2
CHARACTER(LEN=30) :: cact3
CHARACTER(LEN=55) :: cact4

! mppac 

INTEGER, PARAMETER:: inxdim    = 15000

! for Statistics Calcs in DIAGO ; Prints in RDOBS,GETOBS
REAL :: r_stat(model_levels_max,0:8)
REAL :: s_stat(model_levels_max,0:8)
INTEGER :: counta(nobtypmx)
INTEGER :: countb(nobtypmx)
INTEGER :: countc(nobtypmx)

! to pass longitudes and latitudes for edges of local
! box from setcona to RDOBS and HINTCF
REAL :: long_e
REAL :: long_w
REAL :: lat_n,lat_s
REAL :: long_w_model
REAL :: long_e_model

END MODULE ac_control_mod
