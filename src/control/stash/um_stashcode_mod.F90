! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  List of stashcode magic numbers

MODULE um_stashcode_mod

! Description:
!   Stash code definitions used in the RCF, ancillary update and PWS codes
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stash
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

! Section numbers 
INTEGER, PARAMETER :: stashcode_prog_sec          =   0
INTEGER, PARAMETER :: stashcode_bl_sec            =   3
INTEGER, PARAMETER :: stashcode_conv_sec          =   5
INTEGER, PARAMETER :: stashcode_gwd_sec           =   6
INTEGER, PARAMETER :: stashcode_proc_dyn_sec      =  15
INTEGER, PARAMETER :: stashcode_proc_phys_sec     =  16
INTEGER, PARAMETER :: stashcode_aerosol_sec       =  17
INTEGER, PARAMETER :: stashcode_pws_sec           =  20
INTEGER, PARAMETER :: stashcode_lbc_input_sec     =  31 
INTEGER, PARAMETER :: stashcode_lbc_output_sec    =  32 
INTEGER, PARAMETER :: stashcode_tracer_sec        =  33
INTEGER, PARAMETER :: stashcode_ukca_sec          =  34
INTEGER, PARAMETER :: stashcode_tracer_lbc_sec    =  36 
INTEGER, PARAMETER :: stashcode_ukca_lbc_sec      =  37 
INTEGER, PARAMETER :: stashcode_glomap_sec        =  38
INTEGER, PARAMETER :: stashcode_ukca_chem_diag    =  50
INTEGER, PARAMETER :: stashcode_glomap_clim_sec   =  54

! Stashcodes used in the reconfiguration and/or ancillary updating

INTEGER, PARAMETER :: stashcode_u                 =   2
INTEGER, PARAMETER :: stashcode_v                 =   3
INTEGER, PARAMETER :: stashcode_theta             =   4
INTEGER, PARAMETER :: stashcode_orog_x_grad       =   5
INTEGER, PARAMETER :: stashcode_orog_y_grad       =   6
INTEGER, PARAMETER :: stashcode_unfilt_orog       =   7
INTEGER, PARAMETER :: stashcode_soil_moist        =   9

INTEGER, PARAMETER :: stashcode_q                 =  10
INTEGER, PARAMETER :: stashcode_qcf               =  12
INTEGER, PARAMETER :: stashcode_cca               =  13
INTEGER, PARAMETER :: stashcode_ccb               =  14
INTEGER, PARAMETER :: stashcode_cct               =  15
INTEGER, PARAMETER :: stashcode_cc_lwp            =  16
INTEGER, PARAMETER :: stashcode_sil_orog_rough    =  17
INTEGER, PARAMETER :: stashcode_hlf_pk_to_trf_ht  =  18

INTEGER, PARAMETER :: stashcode_soil_temp         =  20
INTEGER, PARAMETER :: stashcode_lcbase            =  21
INTEGER, PARAMETER :: stashcode_mean_canopyw      =  22
INTEGER, PARAMETER :: stashcode_snow_amount       =  23
INTEGER, PARAMETER :: stashcode_mean_snow         =  23
INTEGER, PARAMETER :: stashcode_surftemp          =  24
INTEGER, PARAMETER :: stashcode_tstar             =  24
INTEGER, PARAMETER :: stashcode_bl_depth          =  25
INTEGER, PARAMETER :: stashcode_rough_length      =  26
INTEGER, PARAMETER :: stashcode_z0                =  26
INTEGER, PARAMETER :: stashcode_snow_edge         =  27
INTEGER, PARAMETER :: stashcode_surf_z_curr       =  28
INTEGER, PARAMETER :: stashcode_surf_m_curr       =  29

INTEGER, PARAMETER :: stashcode_lsm               =  30
INTEGER, PARAMETER :: stashcode_icefrac           =  31
INTEGER, PARAMETER :: stashcode_icethick          =  32
INTEGER, PARAMETER :: stashcode_orog              =  33
INTEGER, PARAMETER :: stashcode_orog_var          =  34
INTEGER, PARAMETER :: stashcode_orog_gdxx         =  35
INTEGER, PARAMETER :: stashcode_orog_gdxy         =  36
INTEGER, PARAMETER :: stashcode_orog_gdyy         =  37
INTEGER, PARAMETER :: stashcode_ice_edge_ancil    =  38
INTEGER, PARAMETER :: stashcode_ice_edge_inancil  =  38
INTEGER, PARAMETER :: stashcode_tstar_anom        =  39

INTEGER, PARAMETER :: stashcode_vol_smc_wilt      =  40
INTEGER, PARAMETER :: stashcode_vol_smc_cri       =  41
INTEGER, PARAMETER :: stashcode_vol_smc_sat       =  43
INTEGER, PARAMETER :: stashcode_Ksat              =  44
INTEGER, PARAMETER :: stashcode_thermal_capacity  =  46
INTEGER, PARAMETER :: stashcode_thermal_conduct   =  47
INTEGER, PARAMETER :: stashcode_soil_suction      =  48
INTEGER, PARAMETER :: stashcode_sea_ice_temp      =  49

INTEGER, PARAMETER :: stashcode_veg_frac          =  50
INTEGER, PARAMETER :: stashcode_total_aero_emiss  =  57
INTEGER, PARAMETER :: stashcode_SO2_emiss         =  58
INTEGER, PARAMETER :: stashcode_dimethyl_sul_emiss=  59

INTEGER, PARAMETER :: stashcode_ozone             =  60
INTEGER, PARAMETER :: stashcode_e_trb             =  70
INTEGER, PARAMETER :: stashcode_tsq_trb           =  71
INTEGER, PARAMETER :: stashcode_qsq_trb           =  72
INTEGER, PARAMETER :: stashcode_cov_trb           =  73
INTEGER, PARAMETER :: stashcode_zhpar_shcu        =  74

INTEGER, PARAMETER :: stashcode_cloud_number      =  75
INTEGER, PARAMETER :: stashcode_rain_number       =  76
INTEGER, PARAMETER :: stashcode_rain_3mom         =  77
INTEGER, PARAMETER :: stashcode_ice_number        =  78
INTEGER, PARAMETER :: stashcode_snow_number       =  79
INTEGER, PARAMETER :: stashcode_snow_3mom         =  80
INTEGER, PARAMETER :: stashcode_graup_number      =  81
INTEGER, PARAMETER :: stashcode_graup_3mom        =  82

INTEGER, PARAMETER :: stashcode_activesol_liquid  =  83
INTEGER, PARAMETER :: stashcode_activesol_rain    =  84
INTEGER, PARAMETER :: stashcode_active_insol_ice  =  85
INTEGER, PARAMETER :: stashcode_active_sol_ice    =  86
INTEGER, PARAMETER :: stashcode_active_insol_liq  =  87
INTEGER, PARAMETER :: stashcode_active_sol_num    =  88
INTEGER, PARAMETER :: stashcode_active_insol_num  =  89

INTEGER, PARAMETER :: stashcode_total_aero        =  90
INTEGER, PARAMETER :: stashcode_flash_pot         =  91
INTEGER, PARAMETER :: stashcode_runoff_coast_out  =  93
INTEGER, PARAMETER :: stashcode_snow_on_ice       =  95
INTEGER, PARAMETER :: stashcode_ocnsrf_chlorophyll=  96
INTEGER, PARAMETER :: stashcode_chlorophyll       =  96
INTEGER, PARAMETER :: stashcode_z0m_soil          =  97

INTEGER, PARAMETER :: stashcode_blwvariance       =  99

INTEGER, PARAMETER :: stashcode_so2               = 101
INTEGER, PARAMETER :: stashcode_dms               = 102
INTEGER, PARAMETER :: stashcode_mmr_so4_aitken    = 103
INTEGER, PARAMETER :: stashcode_mmr_so4_accum     = 104
INTEGER, PARAMETER :: stashcode_mmr_so4_diss      = 105
INTEGER, PARAMETER :: stashcode_mmr_nh3           = 107
INTEGER, PARAMETER :: stashcode_mmr_bc_fr         = 108
INTEGER, PARAMETER :: stashcode_mmr_bc_ag         = 109
INTEGER, PARAMETER :: stashcode_mmr_bc_cl         = 110
INTEGER, PARAMETER :: stashcode_mmr_smoke_fr      = 111
INTEGER, PARAMETER :: stashcode_mmr_smoke_ag      = 112
INTEGER, PARAMETER :: stashcode_mmr_smoke_cl      = 113
INTEGER, PARAMETER :: stashcode_mmr_ocff_fr       = 114
INTEGER, PARAMETER :: stashcode_mmr_ocff_ag       = 115
INTEGER, PARAMETER :: stashcode_mmr_ocff_cl       = 116
INTEGER, PARAMETER :: stashcode_mmr_nitr_acc      = 117
INTEGER, PARAMETER :: stashcode_mmr_nitr_diss     = 118

INTEGER, PARAMETER :: stashcode_biom_elev_em_h1   = 119
INTEGER, PARAMETER :: stashcode_biom_elev_em_h2   = 120
INTEGER, PARAMETER :: stashcode_3d_nat_so2_em     = 121
INTEGER, PARAMETER :: stashcode_3d_oh_conc        = 122
INTEGER, PARAMETER :: stashcode_3d_ho2_conc       = 123
INTEGER, PARAMETER :: stashcode_3dh2o2_mixrat     = 124
INTEGER, PARAMETER :: stashcode_3d_ozone_mixrat   = 125
INTEGER, PARAMETER :: stashcode_hi_SO2_emiss_emiss= 126
INTEGER, PARAMETER :: stashcode_ammonia_gas_emiss = 127
INTEGER, PARAMETER :: stashcode_soot_surf         = 128
INTEGER, PARAMETER :: stashcode_soot_hi_lev       = 129

INTEGER, PARAMETER :: stashcode_biom_surf_em      = 130
INTEGER, PARAMETER :: stashcode_biom_elev_em      = 131
INTEGER, PARAMETER :: stashcode_dms_conc_sea      = 132
INTEGER, PARAMETER :: stashcode_dms_conc_sw       = 132
INTEGER, PARAMETER :: stashcode_ocff_surf_emiss   = 134
INTEGER, PARAMETER :: stashcode_ocff_hilev_emiss  = 135

INTEGER, PARAMETER :: stashcode_w                 = 150
INTEGER, PARAMETER :: stashcode_riv_sequence      = 151
INTEGER, PARAMETER :: stashcode_riv_direction     = 152
INTEGER, PARAMETER :: stashcode_riv_storage       = 153
INTEGER, PARAMETER :: stashcode_ice_subl_cat      = 182
INTEGER, PARAMETER :: stashcode_iceberg_calving   = 190
INTEGER, PARAMETER :: stashcode_sstfrz            = 194
INTEGER, PARAMETER :: stashcode_tstar_ice_cat_cpl = 195

INTEGER, PARAMETER :: stashcode_u_compnt_pert     = 202
INTEGER, PARAMETER :: stashcode_v_compnt_pert     = 203
INTEGER, PARAMETER :: stashcode_clapp_hb          = 207
INTEGER, PARAMETER :: stashcode_3d_cca            = 211
INTEGER, PARAMETER :: stashcode_3d_ccw            = 212
INTEGER, PARAMETER :: stashcode_can_conduct       = 213
INTEGER, PARAMETER :: stashcode_unfrozen_soil     = 214
INTEGER, PARAMETER :: stashcode_frozen_soil       = 215
INTEGER, PARAMETER :: stashcode_frac_surf_type    = 216
INTEGER, PARAMETER :: stashcode_lai               = 217
INTEGER, PARAMETER :: stashcode_canopy_height     = 218
INTEGER, PARAMETER :: stashcode_disturb_frac_veg  = 219

INTEGER, PARAMETER :: stashcode_snw_free_alb_bs   = 220
INTEGER, PARAMETER :: stashcode_soil_carbon_cont  = 223
INTEGER, PARAMETER :: stashcode_npp_pft_acc       = 224
INTEGER, PARAMETER :: stashcode_g_lf_pft_acc      = 225
INTEGER, PARAMETER :: stashcode_g_ph_lf_pft_acc   = 226
INTEGER, PARAMETER :: stashcode_rsp_w_pft_acc     = 227
INTEGER, PARAMETER :: stashcode_rsp_s_acc         = 228
INTEGER, PARAMETER :: stashcode_can_water_tile    = 229

INTEGER, PARAMETER :: stashcode_catch_tile        = 230
INTEGER, PARAMETER :: stashcode_rgrain            = 231
INTEGER, PARAMETER :: stashcode_tstar_tile        = 233
INTEGER, PARAMETER :: stashcode_tsurf_elev_surft  = 576
INTEGER, PARAMETER :: stashcode_z0_tile           = 234
INTEGER, PARAMETER :: stashcode_infil_max_tile    = 236
INTEGER, PARAMETER :: stashcode_sw_down_tile      = 237
INTEGER, PARAMETER :: stashcode_sw_down           = 238
INTEGER, PARAMETER :: stashcode_lw_up_diff        = 239

INTEGER, PARAMETER :: stashcode_snow_tile         = 240
INTEGER, PARAMETER :: stashcode_catch_snow        = 241
INTEGER, PARAMETER :: stashcode_snow_grnd         = 242
INTEGER, PARAMETER :: stashcode_surf_sw_alb       = 243
INTEGER, PARAMETER :: stashcode_surf_vis_alb      = 244
INTEGER, PARAMETER :: stashcode_surf_nir_alb      = 245
INTEGER, PARAMETER :: stashcode_z0h_tile          = 246

INTEGER, PARAMETER :: stashcode_CO2_surf_emiss    = 251
INTEGER, PARAMETER :: stashcode_rho               = 253
INTEGER, PARAMETER :: stashcode_qcl               = 254
INTEGER, PARAMETER :: stashcode_exner             = 255
INTEGER, PARAMETER :: stashcode_u_adv             = 256
INTEGER, PARAMETER :: stashcode_v_adv             = 257
INTEGER, PARAMETER :: stashcode_w_adv             = 258
INTEGER, PARAMETER :: stashcode_n_turb_mixlvs     = 259

INTEGER, PARAMETER :: stashcode_lvl_bse_dp_sc     = 260
INTEGER, PARAMETER :: stashcode_lvl_top_dp_sc     = 261
INTEGER, PARAMETER :: stashcode_bl_conv_flag      = 262
INTEGER, PARAMETER :: stashcode_turb_temp         = 263
INTEGER, PARAMETER :: stashcode_turb_humid        = 264
INTEGER, PARAMETER :: stashcode_area_cf           = 265
INTEGER, PARAMETER :: stashcode_bulk_cf           = 266
INTEGER, PARAMETER :: stashcode_liquid_cf         = 267
INTEGER, PARAMETER :: stashcode_frozen_cf         = 268
INTEGER, PARAMETER :: stashcode_sfc_zonal_cur     = 269

INTEGER, PARAMETER :: stashcode_sfc_merid_cur     = 270
INTEGER, PARAMETER :: stashcode_qcf2              = 271
INTEGER, PARAMETER :: stashcode_qrain             = 272
INTEGER, PARAMETER :: stashcode_qgraup            = 273
INTEGER, PARAMETER :: stashcode_top_ind_mean      = 274
INTEGER, PARAMETER :: stashcode_Ti_Mean           = 274
INTEGER, PARAMETER :: stashcode_top_ind_stddev    = 275
INTEGER, PARAMETER :: stashcode_Ti_Sig            = 275
INTEGER, PARAMETER :: stashcode_fexp              = 276
INTEGER, PARAMETER :: stashcode_gamtot            = 277
INTEGER, PARAMETER :: stashcode_zw                = 278
INTEGER, PARAMETER :: stashcode_fsat              = 279

INTEGER, PARAMETER :: stashcode_fwetl             = 280
INTEGER, PARAMETER :: stashcode_sthzw             = 281
INTEGER, PARAMETER :: stashcode_a_fsat            = 282
INTEGER, PARAMETER :: stashcode_c_fsat            = 283
INTEGER, PARAMETER :: stashcode_a_fwet            = 284
INTEGER, PARAMETER :: stashcode_c_fwet            = 285

INTEGER, PARAMETER :: stashcode_disturb_frac_veg_prev = 286
INTEGER, PARAMETER :: stashcode_wood_prod_fast    = 287
INTEGER, PARAMETER :: stashcode_wood_prod_med     = 288
INTEGER, PARAMETER :: stashcode_wood_prod_slow    = 289

INTEGER, PARAMETER :: stashcode_flake_depth       = 291
INTEGER, PARAMETER :: stashcode_flake_fetch       = 292
INTEGER, PARAMETER :: stashcode_flake_t_mean      = 293
INTEGER, PARAMETER :: stashcode_flake_t_mxl       = 294
INTEGER, PARAMETER :: stashcode_flake_t_ice       = 295
INTEGER, PARAMETER :: stashcode_flake_h_mxl       = 296
INTEGER, PARAMETER :: stashcode_flake_h_ice       = 297
INTEGER, PARAMETER :: stashcode_flake_shape       = 298
INTEGER, PARAMETER :: stashcode_flake_g_over_dt   = 299

INTEGER, PARAMETER :: stashcode_user_anc_sing1    = 301
INTEGER, PARAMETER :: stashcode_user_anc_sing20   = 320
INTEGER, PARAMETER :: stashcode_user_anc_mult1    = 321

INTEGER, PARAMETER :: stashcode_user_anc_mult20   = 340
INTEGER, PARAMETER :: stashcode_tppsozone         = 341
INTEGER, PARAMETER :: stashcode_deep_conv_flag    = 342
INTEGER, PARAMETER :: stashcode_past_conv_precip  = 343
INTEGER, PARAMETER :: stashcode_past_conv_depth   = 344
INTEGER, PARAMETER :: stashcode_cca_dp            = 345
INTEGER, PARAMETER :: stashcode_cca_md            = 346
INTEGER, PARAMETER :: stashcode_cca_sh            = 347
INTEGER, PARAMETER :: stashcode_total_precip      = 348

INTEGER, PARAMETER :: stashcode_clim_biogenic_aero= 351
INTEGER, PARAMETER :: stashcode_clim_delta_aero   = 371
INTEGER, PARAMETER :: stashcode_snowdep_grd_tile  = 376
INTEGER, PARAMETER :: stashcode_snowpack_bk_dens  = 377

INTEGER, PARAMETER :: stashcode_nsnow_layrs_tiles = 380
INTEGER, PARAMETER :: stashcode_snow_laythk_tiles = 381
INTEGER, PARAMETER :: stashcode_snow_ice_tile     = 382
INTEGER, PARAMETER :: stashcode_snow_liq_tile     = 383
INTEGER, PARAMETER :: stashcode_snow_T_tile       = 384
INTEGER, PARAMETER :: stashcode_snow_laydns_tiles = 385
INTEGER, PARAMETER :: stashcode_snow_grnsiz_tiles = 386

INTEGER, PARAMETER :: stashcode_p                 = 407
INTEGER, PARAMETER :: stashcode_pstar             = 409
INTEGER, PARAMETER :: stashcode_ice_conc_cat      = 413
INTEGER, PARAMETER :: stashcode_ice_thick_cat     = 414
INTEGER, PARAMETER :: stashcode_ice_temp_cat      = 415
INTEGER, PARAMETER :: stashcode_ice_snow_depth    = 416
INTEGER, PARAMETER :: stashcode_dust_parent_clay  = 418
INTEGER, PARAMETER :: stashcode_dust_parent_silt  = 419

INTEGER, PARAMETER :: stashcode_dust_parent_sand  = 420
INTEGER, PARAMETER :: stashcode_dust_soil_mf1     = 421
INTEGER, PARAMETER :: stashcode_dust_soil_mf2     = 422
INTEGER, PARAMETER :: stashcode_dust_soil_mf3     = 423
INTEGER, PARAMETER :: stashcode_dust_soil_mf4     = 424
INTEGER, PARAMETER :: stashcode_dust_soil_mf5     = 425
INTEGER, PARAMETER :: stashcode_dust_soil_mf6     = 426
INTEGER, PARAMETER :: stashcode_soil_massfrac6    = 426
INTEGER, PARAMETER :: stashcode_pond_frac_cat     = 428
INTEGER, PARAMETER :: stashcode_pond_depth_cat    = 429

INTEGER, PARAMETER :: stashcode_dust1_mmr         = 431
INTEGER, PARAMETER :: stashcode_dust2_mmr         = 432
INTEGER, PARAMETER :: stashcode_dust3_mmr         = 433
INTEGER, PARAMETER :: stashcode_dust4_mmr         = 434
INTEGER, PARAMETER :: stashcode_dust5_mmr         = 435
INTEGER, PARAMETER :: stashcode_dust6_mmr         = 436

INTEGER, PARAMETER :: stashcode_ice_surf_cond_cat = 440
INTEGER, PARAMETER :: stashcode_ice_surf_temp_cat = 441

INTEGER, PARAMETER :: stashcode_soilnitro_dpm     = 442
INTEGER, PARAMETER :: stashcode_soilnitro_rpm     = 443
INTEGER, PARAMETER :: stashcode_soilnitro_bio     = 444
INTEGER, PARAMETER :: stashcode_soilnitro_hum     = 445
INTEGER, PARAMETER :: stashcode_soil_inorgnit     = 446
INTEGER, PARAMETER :: stashcode_nitrogen_deposition = 447

INTEGER, PARAMETER :: stashcode_crop_frac         = 448
INTEGER, PARAMETER :: stashcode_pasture_frac      = 458

INTEGER, PARAMETER :: stashcode_soilcarb_dpm      = 466
INTEGER, PARAMETER :: stashcode_soilcarb_rpm      = 467
INTEGER, PARAMETER :: stashcode_soilcarb_bio      = 468
INTEGER, PARAMETER :: stashcode_soilcarb_hum      = 469

INTEGER, PARAMETER :: stashcode_ozone_tracer      = 480
INTEGER, PARAMETER :: stashcode_o3_prod_loss      = 481
INTEGER, PARAMETER :: stashcode_o3_p_l_vmr        = 482
INTEGER, PARAMETER :: stashcode_o3_vmr            = 483
INTEGER, PARAMETER :: stashcode_o3_p_l_temp       = 484
INTEGER, PARAMETER :: stashcode_o3_temp           = 485
INTEGER, PARAMETER :: stashcode_o3_p_l_colo3      = 486
INTEGER, PARAMETER :: stashcode_o3_colo3          = 487

INTEGER, PARAMETER :: stashcode_dctemp_tile       = 490
INTEGER, PARAMETER :: stashcode_dctemp_ssi        = 491
INTEGER, PARAMETER :: stashcode_tm_trans          = 492
INTEGER, PARAMETER :: stashcode_ddmfx             = 493
INTEGER, PARAMETER :: stashcode_urbhgt            = 494
INTEGER, PARAMETER :: stashcode_urbhwr            = 495
INTEGER, PARAMETER :: stashcode_urbwrr            = 496
INTEGER, PARAMETER :: stashcode_urbdisp           = 497
INTEGER, PARAMETER :: stashcode_urbztm            = 498
INTEGER, PARAMETER :: stashcode_urbalbwl          = 499

INTEGER, PARAMETER :: stashcode_urbalbrd          = 500
INTEGER, PARAMETER :: stashcode_urbemisw          = 501
INTEGER, PARAMETER :: stashcode_urbemisr          = 502
INTEGER, PARAMETER :: stashcode_land_frac         = 505
INTEGER, PARAMETER :: stashcode_tstar_land        = 506
INTEGER, PARAMETER :: stashcode_tstar_sea         = 507
INTEGER, PARAMETER :: stashcode_tstar_sice        = 508
INTEGER, PARAMETER :: stashcode_albedo_sice       = 509
INTEGER, PARAMETER :: stashcode_albedo_land       = 510

INTEGER, PARAMETER :: stashcode_ux_ccp            = 569
INTEGER, PARAMETER :: stashcode_uy_ccp            = 570
INTEGER, PARAMETER :: stashcode_um_ccp            = 571
INTEGER, PARAMETER :: stashcode_g_ccp             = 572
INTEGER, PARAMETER :: stashcode_h_ccp             = 573
INTEGER, PARAMETER :: stashcode_riso_ccp          = 574
INTEGER, PARAMETER :: stashcode_rdir_ccp          = 575

! PV-tracers 
INTEGER, PARAMETER :: stashcode_dPV_rad           = 577
INTEGER, PARAMETER :: stashcode_dPV_sw            = 578
INTEGER, PARAMETER :: stashcode_dPV_lw            = 579
INTEGER, PARAMETER :: stashcode_dPV_mic           = 580
INTEGER, PARAMETER :: stashcode_dPV_gwd           = 581
INTEGER, PARAMETER :: stashcode_dPV_ph1           = 582
INTEGER, PARAMETER :: stashcode_dPV_conv          = 583
INTEGER, PARAMETER :: stashcode_dPV_bl            = 584
INTEGER, PARAMETER :: stashcode_dPV_stph          = 585
INTEGER, PARAMETER :: stashcode_dPV_cld           = 586
INTEGER, PARAMETER :: stashcode_dPV_iau           = 587
INTEGER, PARAMETER :: stashcode_dPV_nud           = 588
INTEGER, PARAMETER :: stashcode_dPV_tot           = 589
INTEGER, PARAMETER :: stashcode_dPV_adv           = 590
INTEGER, PARAMETER :: stashcode_dPV_sol           = 591
INTEGER, PARAMETER :: stashcode_dPV_mass          = 592
INTEGER, PARAMETER :: stashcode_adv_only_PV       = 593

! Diabatic tracers
INTEGER, PARAMETER :: stashcode_dtheta_0          = 600
INTEGER, PARAMETER :: stashcode_dtheta_bl         = 601
INTEGER, PARAMETER :: stashcode_dtheta_bl_mix     = 602
INTEGER, PARAMETER :: stashcode_dtheta_bl_LH      = 603
INTEGER, PARAMETER :: stashcode_dtheta_conv       = 604
INTEGER, PARAMETER :: stashcode_dtheta_mic        = 605
INTEGER, PARAMETER :: stashcode_dtheta_rad        = 606
INTEGER, PARAMETER :: stashcode_dtheta_SW         = 607
INTEGER, PARAMETER :: stashcode_dtheta_LW         = 608
INTEGER, PARAMETER :: stashcode_dtheta_slow       = 609
INTEGER, PARAMETER :: stashcode_dtheta_cld        = 610

! Stochastic Physics
INTEGER, PARAMETER :: stashcode_bl_pert_rand_fld  = 595
INTEGER, PARAMETER :: stashcode_bl_pert_flag      = 596

!----------------------------------------------------------
! Section 16 fields that may be reconfigured for VAR
!----------------------------------------------------------
INTEGER, PARAMETER :: stashcode_t                 = 16004
INTEGER, PARAMETER :: stashcode_qc                = 16206
INTEGER, PARAMETER :: stashcode_qT                = 16207

!----------------------------------------------------------
! UKCA stashcodes - section 34
!----------------------------------------------------------
INTEGER, PARAMETER :: stashcode_NO2               = 34004
INTEGER, PARAMETER :: stashcode_CH4               = 34009
INTEGER, PARAMETER :: stashcode_CO                = 34010
INTEGER, PARAMETER :: stashcode_HCHO              = 34011
INTEGER, PARAMETER :: stashcode_O3                = 34001
INTEGER, PARAMETER :: stashcode_NO                = 34002
INTEGER, PARAMETER :: stashcode_HNO3              = 34007
INTEGER, PARAMETER :: stashcode_PAN               = 34017
INTEGER, PARAMETER :: stashcode_C2H6              = 34014
INTEGER, PARAMETER :: stashcode_C3H8              = 34018
! Aerosols
INTEGER, PARAMETER :: stashcode_nucsol_no         = 34101
INTEGER, PARAMETER :: stashcode_nucsol_so4        = 34102
INTEGER, PARAMETER :: stashcode_Aitsol_no         = 34103
INTEGER, PARAMETER :: stashcode_Aitsol_so4        = 34104
INTEGER, PARAMETER :: stashcode_Aitsol_bc         = 34105
INTEGER, PARAMETER :: stashcode_Aitsol_oc         = 34106
INTEGER, PARAMETER :: stashcode_accsol_no         = 34107
INTEGER, PARAMETER :: stashcode_accsol_so4        = 34108
INTEGER, PARAMETER :: stashcode_accsol_bc         = 34109
INTEGER, PARAMETER :: stashcode_accsol_oc         = 34110
INTEGER, PARAMETER :: stashcode_accsol_ss         = 34111
INTEGER, PARAMETER :: stashcode_accsol_du         = 34112
INTEGER, PARAMETER :: stashcode_corsol_no         = 34113
INTEGER, PARAMETER :: stashcode_corsol_so4        = 34114
INTEGER, PARAMETER :: stashcode_corsol_bc         = 34115
INTEGER, PARAMETER :: stashcode_corsol_oc         = 34116
INTEGER, PARAMETER :: stashcode_corsol_ss         = 34117
INTEGER, PARAMETER :: stashcode_corsol_du         = 34118
INTEGER, PARAMETER :: stashcode_Aitinsol_no       = 34119
INTEGER, PARAMETER :: stashcode_Aitinsol_bc       = 34120
INTEGER, PARAMETER :: stashcode_Aitinsol_oc       = 34121
INTEGER, PARAMETER :: stashcode_accinsol_no       = 34122
INTEGER, PARAMETER :: stashcode_accinsol_du       = 34123
INTEGER, PARAMETER :: stashcode_corinsol_no       = 34124
INTEGER, PARAMETER :: stashcode_corinsol_du       = 34125
INTEGER, PARAMETER :: stashcode_nucsol_oc         = 34126
INTEGER, PARAMETER :: stashcode_nucsol_so         = 34128
INTEGER, PARAMETER :: stashcode_Aitsol_so         = 34129
INTEGER, PARAMETER :: stashcode_accsol_so         = 34130
INTEGER, PARAMETER :: stashcode_corsol_so         = 34131
INTEGER, PARAMETER :: stashcode_accsol_nh4        = 34134
INTEGER, PARAMETER :: stashcode_corsol_nh4        = 34135

! RADAER prognostics
INTEGER, PARAMETER :: stashcode_dryd_ait_sol      = 34921
INTEGER, PARAMETER :: stashcode_dryd_acc_sol      = 34922
INTEGER, PARAMETER :: stashcode_dryd_cor_sol      = 34923
INTEGER, PARAMETER :: stashcode_dryd_ait_insol    = 34924
INTEGER, PARAMETER :: stashcode_dryd_acc_insol    = 34925
INTEGER, PARAMETER :: stashcode_dryd_cor_insol    = 34926
INTEGER, PARAMETER :: stashcode_wetd_ait_sol      = 34927
INTEGER, PARAMETER :: stashcode_wetd_acc_sol      = 34928
INTEGER, PARAMETER :: stashcode_wetd_cor_sol      = 34929
INTEGER, PARAMETER :: stashcode_rho_ait_sol       = 34930
INTEGER, PARAMETER :: stashcode_rho_acc_sol       = 34931
INTEGER, PARAMETER :: stashcode_rho_cor_sol       = 34932
INTEGER, PARAMETER :: stashcode_rho_ait_insol     = 34933
INTEGER, PARAMETER :: stashcode_rho_acc_insol     = 34934
INTEGER, PARAMETER :: stashcode_rho_cor_insol     = 34935
INTEGER, PARAMETER :: stashcode_pvol_ait_su_sol   = 34936
INTEGER, PARAMETER :: stashcode_pvol_ait_bc_sol   = 34937
INTEGER, PARAMETER :: stashcode_pvol_ait_oc_sol   = 34938
INTEGER, PARAMETER :: stashcode_pvol_ait_no3_sol  = 34939
INTEGER, PARAMETER :: stashcode_pvol_ait_so_sol   = 34940
INTEGER, PARAMETER :: stashcode_pvol_ait_h2o_sol  = 34941
INTEGER, PARAMETER :: stashcode_pvol_acc_su_sol   = 34942
INTEGER, PARAMETER :: stashcode_pvol_acc_bc_sol   = 34943
INTEGER, PARAMETER :: stashcode_pvol_acc_oc_sol   = 34944
INTEGER, PARAMETER :: stashcode_pvol_acc_ss_sol   = 34945
INTEGER, PARAMETER :: stashcode_pvol_acc_no3_sol  = 34946
INTEGER, PARAMETER :: stashcode_pvol_acc_du_sol   = 34947
INTEGER, PARAMETER :: stashcode_pvol_acc_so_sol   = 34948
INTEGER, PARAMETER :: stashcode_pvol_acc_h2o_sol  = 34951
INTEGER, PARAMETER :: stashcode_pvol_cor_su_sol   = 34952
INTEGER, PARAMETER :: stashcode_pvol_cor_bc_sol   = 34953
INTEGER, PARAMETER :: stashcode_pvol_cor_oc_sol   = 34954
INTEGER, PARAMETER :: stashcode_pvol_cor_ss_sol   = 34955
INTEGER, PARAMETER :: stashcode_pvol_cor_no3_sol  = 34956
INTEGER, PARAMETER :: stashcode_pvol_cor_du_sol   = 34957
INTEGER, PARAMETER :: stashcode_pvol_cor_so_sol   = 34958
INTEGER, PARAMETER :: stashcode_pvol_cor_h2o_sol  = 34961
INTEGER, PARAMETER :: stashcode_pvol_ait_bc_insol = 34962
INTEGER, PARAMETER :: stashcode_pvol_ait_oc_insol = 34963
INTEGER, PARAMETER :: stashcode_pvol_acc_du_insol = 34964
INTEGER, PARAMETER :: stashcode_pvol_cor_du_insol = 34965

! Non transported prognostics
INTEGER, PARAMETER :: stashcode_surfarea      = 34966
INTEGER, PARAMETER :: stashcode_cdnc3         = 34967
INTEGER, PARAMETER :: stashcode_cdnc          = 34968
INTEGER, PARAMETER :: stashcode_ho2s          = 34969
INTEGER, PARAMETER :: stashcode_ohs           = 34970
INTEGER, PARAMETER :: stashcode_o1ds          = 34971
INTEGER, PARAMETER :: stashcode_o3ps          = 34972
INTEGER, PARAMETER :: stashcode_het_ho2       = 34973
INTEGER, PARAMETER :: stashcode_het_n2o5      = 34974
INTEGER, PARAMETER :: stashcode_tolp1         = 34975
INTEGER, PARAMETER :: stashcode_hoipo2        = 34976
INTEGER, PARAMETER :: stashcode_homvko2       = 34977
INTEGER, PARAMETER :: stashcode_memald1       = 34978
INTEGER, PARAMETER :: stashcode_oxyl1         = 34979
INTEGER, PARAMETER :: stashcode_hoc3h6o2      = 34980
INTEGER, PARAMETER :: stashcode_hoc2h4o2      = 34981
INTEGER, PARAMETER :: stashcode_meko2         = 34982
INTEGER, PARAMETER :: stashcode_mecoch2oo     = 34983
INTEGER, PARAMETER :: stashcode_mecoc2oo      = 34984
INTEGER, PARAMETER :: stashcode_etco3         = 34985
INTEGER, PARAMETER :: stashcode_iproo         = 34986
INTEGER, PARAMETER :: stashcode_sbuoo         = 34987
INTEGER, PARAMETER :: stashcode_nproo         = 34988
INTEGER, PARAMETER :: stashcode_meco3         = 34989
INTEGER, PARAMETER :: stashcode_etoo          = 34990
INTEGER, PARAMETER :: stashcode_meoo          = 34991
INTEGER, PARAMETER :: stashcode_hcl_unlmp     = 34992
INTEGER, PARAMETER :: stashcode_ho2_ntp       = 34993
INTEGER, PARAMETER :: stashcode_bro_unlmp     = 34994
INTEGER, PARAMETER :: stashcode_oh_ntp        = 34995
INTEGER, PARAMETER :: stashcode_no2_unlmp     = 34996
INTEGER, PARAMETER :: stashcode_o1d_ntp       = 34997
INTEGER, PARAMETER :: stashcode_o3p_ntp       = 34998

!----------------------------------------------------------
! LBC input stashcodes 
!----------------------------------------------------------

INTEGER, PARAMETER :: stashcode_lbc_orog           = 31001
INTEGER, PARAMETER :: stashcode_lbc_u              = 31002 
INTEGER, PARAMETER :: stashcode_lbc_v              = 31003
INTEGER, PARAMETER :: stashcode_lbc_w              = 31004
INTEGER, PARAMETER :: stashcode_lbc_density        = 31005
INTEGER, PARAMETER :: stashcode_lbc_theta          = 31006
INTEGER, PARAMETER :: stashcode_lbc_q              = 31007
INTEGER, PARAMETER :: stashcode_lbc_qcl            = 31008
INTEGER, PARAMETER :: stashcode_lbc_qcf            = 31009
INTEGER, PARAMETER :: stashcode_lbc_exner          = 31010

INTEGER, PARAMETER :: stashcode_lbc_u_adv          = 31011
INTEGER, PARAMETER :: stashcode_lbc_v_adv          = 31012
INTEGER, PARAMETER :: stashcode_lbc_w_adv          = 31013
INTEGER, PARAMETER :: stashcode_lbc_qcf2           = 31014
INTEGER, PARAMETER :: stashcode_lbc_qrain          = 31015
INTEGER, PARAMETER :: stashcode_lbc_qgraup         = 31016
INTEGER, PARAMETER :: stashcode_lbc_cf_bulk        = 31017
INTEGER, PARAMETER :: stashcode_lbc_cf_liquid      = 31018
INTEGER, PARAMETER :: stashcode_lbc_cf_frozen      = 31019
INTEGER, PARAMETER :: stashcode_lbc_murk           = 31020

INTEGER, PARAMETER :: stashcode_lbc_dust1_mmr      = 31023
INTEGER, PARAMETER :: stashcode_lbc_dust2_mmr      = 31024
INTEGER, PARAMETER :: stashcode_lbc_dust3_mmr      = 31025
INTEGER, PARAMETER :: stashcode_lbc_dust4_mmr      = 31026
INTEGER, PARAMETER :: stashcode_lbc_dust5_mmr      = 31027
INTEGER, PARAMETER :: stashcode_lbc_dust6_mmr      = 31028
INTEGER, PARAMETER :: stashcode_lbc_so2            = 31029
INTEGER, PARAMETER :: stashcode_lbc_dms            = 31030

INTEGER, PARAMETER :: stashcode_lbc_so4_aitken     = 31031
INTEGER, PARAMETER :: stashcode_lbc_so4_accu       = 31032
INTEGER, PARAMETER :: stashcode_lbc_so4_diss       = 31033
INTEGER, PARAMETER :: stashcode_lbc_nh3            = 31035
INTEGER, PARAMETER :: stashcode_lbc_soot_new       = 31036
INTEGER, PARAMETER :: stashcode_lbc_soot_agd       = 31037
INTEGER, PARAMETER :: stashcode_lbc_soot_cld       = 31038
INTEGER, PARAMETER :: stashcode_lbc_bmass_new      = 31039
INTEGER, PARAMETER :: stashcode_lbc_bmass_agd      = 31040

INTEGER, PARAMETER :: stashcode_lbc_bmass_cld      = 31041
INTEGER, PARAMETER :: stashcode_lbc_ocff_new       = 31042
INTEGER, PARAMETER :: stashcode_lbc_ocff_agd       = 31043
INTEGER, PARAMETER :: stashcode_lbc_ocff_cld       = 31044
INTEGER, PARAMETER :: stashcode_lbc_nitr_acc       = 31045
INTEGER, PARAMETER :: stashcode_lbc_nitr_diss      = 31046

!----------------------------------------------------------
! LBC input stashcodes - Free Tracers
!----------------------------------------------------------

INTEGER, PARAMETER :: stashcode_lbc_free_tracer_1   = 36001
INTEGER, PARAMETER :: stashcode_lbc_free_tracer_150 = 36150

!----------------------------------------------------------
! LBC input stashcodes - UKCA Tracers
!----------------------------------------------------------

INTEGER, PARAMETER :: stashcode_lbc_ukca_1         = 37001
INTEGER, PARAMETER :: stashcode_lbc_ukca_150       = 37150

!----------------------------------------------------------
! UKCA/GLOMAP diagnostic stashcodes - section 38
!----------------------------------------------------------

! CMIP6 diagnostics for component and number densities
INTEGER, PARAMETER :: stashcode_so4_nuc_sol   = 38485
INTEGER, PARAMETER :: stashcode_so4_ait_sol   = 38486
INTEGER, PARAMETER :: stashcode_so4_acc_sol   = 38487
INTEGER, PARAMETER :: stashcode_so4_cor_sol   = 38488
INTEGER, PARAMETER :: stashcode_bc_ait_sol    = 38489
INTEGER, PARAMETER :: stashcode_bc_acc_sol    = 38490
INTEGER, PARAMETER :: stashcode_bc_cor_sol    = 38491
INTEGER, PARAMETER :: stashcode_bc_ait_insol  = 38492
INTEGER, PARAMETER :: stashcode_oc_nuc_sol    = 38493
INTEGER, PARAMETER :: stashcode_oc_ait_sol    = 38494
INTEGER, PARAMETER :: stashcode_oc_acc_sol    = 38495
INTEGER, PARAMETER :: stashcode_oc_cor_sol    = 38496
INTEGER, PARAMETER :: stashcode_oc_ait_insol  = 38497
INTEGER, PARAMETER :: stashcode_ss_acc_sol    = 38498
INTEGER, PARAMETER :: stashcode_ss_cor_sol    = 38499
INTEGER, PARAMETER :: stashcode_du_acc_sol    = 38500
INTEGER, PARAMETER :: stashcode_du_cor_sol    = 38501
INTEGER, PARAMETER :: stashcode_du_acc_insol  = 38502
INTEGER, PARAMETER :: stashcode_du_cor_insol  = 38503
INTEGER, PARAMETER :: stashcode_n_nuc_sol     = 38504
INTEGER, PARAMETER :: stashcode_n_ait_sol     = 38505
INTEGER, PARAMETER :: stashcode_n_acc_sol     = 38506
INTEGER, PARAMETER :: stashcode_n_cor_sol     = 38507
INTEGER, PARAMETER :: stashcode_n_ait_insol   = 38508
INTEGER, PARAMETER :: stashcode_n_acc_insol   = 38509
INTEGER, PARAMETER :: stashcode_n_cor_insol   = 38510
INTEGER, PARAMETER :: stashcode_h2o_nuc_sol   = 38511
INTEGER, PARAMETER :: stashcode_h2o_ait_sol   = 38512
INTEGER, PARAMETER :: stashcode_h2o_acc_sol   = 38513
INTEGER, PARAMETER :: stashcode_h2o_cor_sol   = 38514
INTEGER, PARAMETER :: stashcode_h2o_total     = 38515
INTEGER, PARAMETER :: stashcode_so4_nuc_sol_load  = 38516
INTEGER, PARAMETER :: stashcode_so4_ait_sol_load  = 38517
INTEGER, PARAMETER :: stashcode_so4_acc_sol_load  = 38518
INTEGER, PARAMETER :: stashcode_so4_cor_sol_load  = 38519
INTEGER, PARAMETER :: stashcode_so4_total_load    = 38520
INTEGER, PARAMETER :: stashcode_bc_ait_sol_load   = 38521
INTEGER, PARAMETER :: stashcode_bc_acc_sol_load   = 38522
INTEGER, PARAMETER :: stashcode_bc_cor_sol_load   = 38523
INTEGER, PARAMETER :: stashcode_bc_ait_insol_load = 38524
INTEGER, PARAMETER :: stashcode_bc_total_load     = 38525
INTEGER, PARAMETER :: stashcode_oc_nuc_sol_load   = 38526
INTEGER, PARAMETER :: stashcode_oc_ait_sol_load   = 38527
INTEGER, PARAMETER :: stashcode_oc_acc_sol_load   = 38528
INTEGER, PARAMETER :: stashcode_oc_cor_sol_load   = 38529
INTEGER, PARAMETER :: stashcode_oc_ait_insol_load = 38530
INTEGER, PARAMETER :: stashcode_oc_total_load     = 38531
INTEGER, PARAMETER :: stashcode_du_acc_sol_load   = 38532
INTEGER, PARAMETER :: stashcode_du_cor_sol_load   = 38533
INTEGER, PARAMETER :: stashcode_du_acc_insol_load = 38534
INTEGER, PARAMETER :: stashcode_du_cor_insol_load = 38535
INTEGER, PARAMETER :: stashcode_du_total_load     = 38536
INTEGER, PARAMETER :: stashcode_ss_acc_sol_load   = 38537
INTEGER, PARAMETER :: stashcode_ss_cor_sol_load   = 38538
INTEGER, PARAMETER :: stashcode_ss_total_load     = 38539
INTEGER, PARAMETER :: stashcode_h2o_nuc_sol_load  = 38540
INTEGER, PARAMETER :: stashcode_h2o_ait_sol_load  = 38541
INTEGER, PARAMETER :: stashcode_h2o_acc_sol_load  = 38542
INTEGER, PARAMETER :: stashcode_h2o_cor_sol_load  = 38543
INTEGER, PARAMETER :: stashcode_h2o_total_load    = 38544
INTEGER, PARAMETER :: stashcode_h2o_mmr           = 38545
INTEGER, PARAMETER :: stashcode_so4_nuc_sol_ps    = 38900
INTEGER, PARAMETER :: stashcode_oc_nuc_sol_ps     = 38901
INTEGER, PARAMETER :: stashcode_so_nuc_sol_ps     = 38902
INTEGER, PARAMETER :: stashcode_so4_ait_sol_ps    = 38905
INTEGER, PARAMETER :: stashcode_bc_ait_sol_ps     = 38906
INTEGER, PARAMETER :: stashcode_oc_ait_sol_ps     = 38907
INTEGER, PARAMETER :: stashcode_so_ait_sol_ps     = 38908
INTEGER, PARAMETER :: stashcode_so4_acc_sol_ps    = 38911
INTEGER, PARAMETER :: stashcode_bc_acc_sol_ps     = 38912
INTEGER, PARAMETER :: stashcode_oc_acc_sol_ps     = 38913
INTEGER, PARAMETER :: stashcode_ss_acc_sol_ps     = 38914
INTEGER, PARAMETER :: stashcode_du_acc_sol_ps     = 38916
INTEGER, PARAMETER :: stashcode_so_acc_sol_ps     = 38917
INTEGER, PARAMETER :: stashcode_no3_acc_sol_ps    = 38918
INTEGER, PARAMETER :: stashcode_nh4_acc_sol_ps    = 38919
INTEGER, PARAMETER :: stashcode_so4_cor_sol_ps    = 38920
INTEGER, PARAMETER :: stashcode_bc_cor_sol_ps     = 38921
INTEGER, PARAMETER :: stashcode_oc_cor_sol_ps     = 38922
INTEGER, PARAMETER :: stashcode_ss_cor_sol_ps     = 38923
INTEGER, PARAMETER :: stashcode_du_cor_sol_ps     = 38925
INTEGER, PARAMETER :: stashcode_so_cor_sol_ps     = 38926
INTEGER, PARAMETER :: stashcode_no3_cor_sol_ps    = 38927
INTEGER, PARAMETER :: stashcode_nh4_cor_sol_ps    = 38928
INTEGER, PARAMETER :: stashcode_bc_ait_insol_ps   = 38929
INTEGER, PARAMETER :: stashcode_oc_ait_insol_ps   = 38930
INTEGER, PARAMETER :: stashcode_du_acc_insol_ps   = 38931
INTEGER, PARAMETER :: stashcode_du_cor_insol_ps   = 38932

!PM10 and PM2.5 diagnostics
INTEGER, PARAMETER :: stashcode_pm10_wet  = 38560
INTEGER, PARAMETER :: stashcode_pm2p5_wet = 38561
INTEGER, PARAMETER :: stashcode_pm10_dry  = 38562
INTEGER, PARAMETER :: stashcode_pm2p5_dry = 38563
INTEGER, PARAMETER :: stashcode_pm10_so4  = 38564
INTEGER, PARAMETER :: stashcode_pm2p5_so4 = 38565
INTEGER, PARAMETER :: stashcode_pm10_bc   = 38566
INTEGER, PARAMETER :: stashcode_pm2p5_bc  = 38567
INTEGER, PARAMETER :: stashcode_pm10_oc   = 38568
INTEGER, PARAMETER :: stashcode_pm2p5_oc  = 38569
INTEGER, PARAMETER :: stashcode_pm10_ss   = 38570
INTEGER, PARAMETER :: stashcode_pm2p5_ss  = 38571
INTEGER, PARAMETER :: stashcode_pm10_du   = 38572
INTEGER, PARAMETER :: stashcode_pm2p5_du  = 38573
!First and last PM diagnostics
INTEGER, PARAMETER :: stashcode_pm_first  = 38560
INTEGER, PARAMETER :: stashcode_pm_last   = 38573

!----------------------------------------------------------
! UKCA chemical diagnostics - section 50
!----------------------------------------------------------

INTEGER, PARAMETER :: stashcode_ukca_nat          = 50218
INTEGER, PARAMETER :: stashcode_ukca_trop_ch4     = 50220
INTEGER, PARAMETER :: stashcode_ukca_trop_o3      = 50221
INTEGER, PARAMETER :: stashcode_ukca_trop_oh      = 50222
INTEGER, PARAMETER :: stashcode_ukca_strat_ch4    = 50223
INTEGER, PARAMETER :: stashcode_ukca_strt_ch4_lss = 50226
INTEGER, PARAMETER :: stashcode_ukca_jo1d         = 50228
INTEGER, PARAMETER :: stashcode_ukca_jn2o         = 50229
INTEGER, PARAMETER :: stashcode_ukca_atmos_ch4    = 50231
INTEGER, PARAMETER :: stashcode_ukca_atmos_co     = 50232
INTEGER, PARAMETER :: stashcode_ukca_atmos_n2o    = 50233
INTEGER, PARAMETER :: stashcode_ukca_atmos_cfc12  = 50234
INTEGER, PARAMETER :: stashcode_ukca_atmos_cfc11  = 50235
INTEGER, PARAMETER :: stashcode_ukca_atmos_ch3br  = 50236
INTEGER, PARAMETER :: stashcode_ukca_atmos_h2     = 50237
INTEGER, PARAMETER :: stashcode_ukca_h2o_incr     = 50240
INTEGER, PARAMETER :: stashcode_ukca_jo2          = 50245
INTEGER, PARAMETER :: stashcode_ukca_jo3p         = 50246
INTEGER, PARAMETER :: stashcode_ukca_so4_sad      = 50256

!----------------------------------------------------------
! GLOMAP_CLIM stashcodes - section 54
!----------------------------------------------------------

! Number density and component mass mixing ratios
INTEGER, PARAMETER :: stashcode_gc_nd_nuc_sol     = 54101
INTEGER, PARAMETER :: stashcode_gc_nuc_sol_su     = 54102
INTEGER, PARAMETER :: stashcode_gc_nd_ait_sol     = 54103
INTEGER, PARAMETER :: stashcode_gc_ait_sol_su     = 54104
INTEGER, PARAMETER :: stashcode_gc_ait_sol_bc     = 54105
INTEGER, PARAMETER :: stashcode_gc_ait_sol_oc     = 54106
INTEGER, PARAMETER :: stashcode_gc_nd_acc_sol     = 54107
INTEGER, PARAMETER :: stashcode_gc_acc_sol_su     = 54108
INTEGER, PARAMETER :: stashcode_gc_acc_sol_bc     = 54109
INTEGER, PARAMETER :: stashcode_gc_acc_sol_oc     = 54110
INTEGER, PARAMETER :: stashcode_gc_acc_sol_ss     = 54111
INTEGER, PARAMETER :: stashcode_gc_acc_sol_du     = 54112
INTEGER, PARAMETER :: stashcode_gc_nd_cor_sol     = 54113
INTEGER, PARAMETER :: stashcode_gc_cor_sol_su     = 54114
INTEGER, PARAMETER :: stashcode_gc_cor_sol_bc     = 54115
INTEGER, PARAMETER :: stashcode_gc_cor_sol_oc     = 54116
INTEGER, PARAMETER :: stashcode_gc_cor_sol_ss     = 54117
INTEGER, PARAMETER :: stashcode_gc_cor_sol_du     = 54118
INTEGER, PARAMETER :: stashcode_gc_nd_ait_ins     = 54119
INTEGER, PARAMETER :: stashcode_gc_ait_ins_bc     = 54120
INTEGER, PARAMETER :: stashcode_gc_ait_ins_oc     = 54121
INTEGER, PARAMETER :: stashcode_gc_nd_acc_ins     = 54122
INTEGER, PARAMETER :: stashcode_gc_acc_ins_du     = 54123
INTEGER, PARAMETER :: stashcode_gc_nd_cor_ins     = 54124
INTEGER, PARAMETER :: stashcode_gc_cor_ins_du     = 54125
INTEGER, PARAMETER :: stashcode_gc_nuc_sol_oc     = 54126
INTEGER, PARAMETER :: stashcode_gc_ait_sol_ss     = 54127
INTEGER, PARAMETER :: stashcode_gc_nuc_sol_so     = 54128
INTEGER, PARAMETER :: stashcode_gc_ait_sol_so     = 54129
INTEGER, PARAMETER :: stashcode_gc_acc_sol_so     = 54130
INTEGER, PARAMETER :: stashcode_gc_cor_sol_so     = 54131
INTEGER, PARAMETER :: stashcode_gc_acc_sol_nh4    = 54134
INTEGER, PARAMETER :: stashcode_gc_cor_sol_nh4    = 54135
INTEGER, PARAMETER :: stashcode_gc_acc_sol_no3    = 54138
INTEGER, PARAMETER :: stashcode_gc_cor_sol_no3    = 54139

! Fields required by RADAER
INTEGER, PARAMETER :: stashcode_gc_dryd_ait_sol      = 54921
INTEGER, PARAMETER :: stashcode_gc_dryd_acc_sol      = 54922
INTEGER, PARAMETER :: stashcode_gc_dryd_cor_sol      = 54923
INTEGER, PARAMETER :: stashcode_gc_dryd_ait_ins      = 54924
INTEGER, PARAMETER :: stashcode_gc_dryd_acc_ins      = 54925
INTEGER, PARAMETER :: stashcode_gc_dryd_cor_ins      = 54926
INTEGER, PARAMETER :: stashcode_gc_wetd_ait_sol      = 54927
INTEGER, PARAMETER :: stashcode_gc_wetd_acc_sol      = 54928
INTEGER, PARAMETER :: stashcode_gc_wetd_cor_sol      = 54929
INTEGER, PARAMETER :: stashcode_gc_rho_ait_sol       = 54930
INTEGER, PARAMETER :: stashcode_gc_rho_acc_sol       = 54931
INTEGER, PARAMETER :: stashcode_gc_rho_cor_sol       = 54932
INTEGER, PARAMETER :: stashcode_gc_rho_ait_ins       = 54933
INTEGER, PARAMETER :: stashcode_gc_rho_acc_ins       = 54934
INTEGER, PARAMETER :: stashcode_gc_rho_cor_ins       = 54935
INTEGER, PARAMETER :: stashcode_gc_pvol_ait_su_sol   = 54936
INTEGER, PARAMETER :: stashcode_gc_pvol_ait_bc_sol   = 54937
INTEGER, PARAMETER :: stashcode_gc_pvol_ait_oc_sol   = 54938
INTEGER, PARAMETER :: stashcode_gc_pvol_ait_no_sol   = 54939
INTEGER, PARAMETER :: stashcode_gc_pvol_ait_so_sol   = 54940
INTEGER, PARAMETER :: stashcode_gc_pvol_ait_h2o_sol  = 54941
INTEGER, PARAMETER :: stashcode_gc_pvol_acc_su_sol   = 54942
INTEGER, PARAMETER :: stashcode_gc_pvol_acc_bc_sol   = 54943
INTEGER, PARAMETER :: stashcode_gc_pvol_acc_oc_sol   = 54944
INTEGER, PARAMETER :: stashcode_gc_pvol_acc_ss_sol   = 54945
INTEGER, PARAMETER :: stashcode_gc_pvol_acc_no3_sol  = 54946
INTEGER, PARAMETER :: stashcode_gc_pvol_acc_du_sol   = 54947
INTEGER, PARAMETER :: stashcode_gc_pvol_acc_so_sol   = 54948
INTEGER, PARAMETER :: stashcode_gc_pvol_acc_h2o_sol  = 54951
INTEGER, PARAMETER :: stashcode_gc_pvol_cor_su_sol   = 54952
INTEGER, PARAMETER :: stashcode_gc_pvol_cor_bc_sol   = 54953
INTEGER, PARAMETER :: stashcode_gc_pvol_cor_oc_sol   = 54954
INTEGER, PARAMETER :: stashcode_gc_pvol_cor_ss_sol   = 54955
INTEGER, PARAMETER :: stashcode_gc_pvol_cor_no3_sol  = 54956
INTEGER, PARAMETER :: stashcode_gc_pvol_cor_du_sol   = 54957
INTEGER, PARAMETER :: stashcode_gc_pvol_cor_so_sol   = 54958
INTEGER, PARAMETER :: stashcode_gc_pvol_cor_h2o_sol  = 54961
INTEGER, PARAMETER :: stashcode_gc_pvol_ait_bc_ins   = 54962
INTEGER, PARAMETER :: stashcode_gc_pvol_ait_oc_ins   = 54963
INTEGER, PARAMETER :: stashcode_gc_pvol_acc_du_ins   = 54964
INTEGER, PARAMETER :: stashcode_gc_pvol_cor_du_ins   = 54965

! Fields required by ACTIVATE

INTEGER, PARAMETER :: stashcode_gc_cf_liquid         = 54476
INTEGER, PARAMETER :: stashcode_gc_cdncwt            = 54477
INTEGER, PARAMETER :: stashcode_gc_cdnc3             = 54967
INTEGER, PARAMETER :: stashcode_gc_cdnc              = 54968

!----------------------------------------------------------
! ENDGame stashcodes
!----------------------------------------------------------

INTEGER, PARAMETER :: stashcode_etadot            = 387
INTEGER, PARAMETER :: stashcode_thetavd           = 388
INTEGER, PARAMETER :: stashcode_dry_rho           = 389
INTEGER, PARAMETER :: stashcode_exner_surf        = 398
INTEGER, PARAMETER :: stashcode_psiw_surf         = 390
INTEGER, PARAMETER :: stashcode_psiw_lid          = 397
INTEGER, PARAMETER :: stashcode_mv                = 391
INTEGER, PARAMETER :: stashcode_mcl               = 392
INTEGER, PARAMETER :: stashcode_mcf               = 393
INTEGER, PARAMETER :: stashcode_mr                = 394
INTEGER, PARAMETER :: stashcode_mgr               = 395
INTEGER, PARAMETER :: stashcode_mcf2              = 396

! ---------------------------------------------------------
! PWS stashcodes - section 20
! ---------------------------------------------------------
INTEGER, PARAMETER :: stashcode_pws_thickness500   =  1
INTEGER, PARAMETER :: stashcode_pws_thickness850   =  2
INTEGER, PARAMETER :: stashcode_pws_windspeed10m   =  3
INTEGER, PARAMETER :: stashcode_pws_windspeedplev  =  4
INTEGER, PARAMETER :: stashcode_pws_divergence     =  5
INTEGER, PARAMETER :: stashcode_pws_rel_vorticity  =  6
INTEGER, PARAMETER :: stashcode_pws_mtn_wave_turb  =  7
INTEGER, PARAMETER :: stashcode_pws_conv_cld_dep   =  12
INTEGER, PARAMETER :: stashcode_pws_precip_sym     =  14
INTEGER, PARAMETER :: stashcode_pws_cat_turb       =  16
INTEGER, PARAMETER :: stashcode_pws_max_cat        =  17
INTEGER, PARAMETER :: stashcode_pws_max_cat_press  =  18
INTEGER, PARAMETER :: stashcode_pws_max_wind_ub    =  20
INTEGER, PARAMETER :: stashcode_pws_max_wind_vb    =  21
INTEGER, PARAMETER :: stashcode_pws_max_wind_pb    =  22
INTEGER, PARAMETER :: stashcode_pws_max_wind_icao  =  23
INTEGER, PARAMETER :: stashcode_pws_snow_prob      =  28
INTEGER, PARAMETER :: stashcode_pws_contrail_bot   =  29
INTEGER, PARAMETER :: stashcode_pws_contrail_top   =  30
INTEGER, PARAMETER :: stashcode_pws_thermal_advec  =  31
INTEGER, PARAMETER :: stashcode_pws_freezing_ht    =  33
INTEGER, PARAMETER :: stashcode_pws_freezing_press =  34
INTEGER, PARAMETER :: stashcode_pws_freezing_icao  =  35
INTEGER, PARAMETER :: stashcode_pws_isotherm_ms20_ht    =  36
INTEGER, PARAMETER :: stashcode_pws_isotherm_ms20_press =  37
INTEGER, PARAMETER :: stashcode_pws_isotherm_ms20_icao  =  38
INTEGER, PARAMETER :: stashcode_pws_conv_icao_base =  39
INTEGER, PARAMETER :: stashcode_pws_conv_icao_top  =  40
INTEGER, PARAMETER :: stashcode_pws_max_wind_base  =  41
INTEGER, PARAMETER :: stashcode_pws_max_wind_top   =  42
INTEGER, PARAMETER :: stashcode_pws_icing_pot_diag =  43
INTEGER, PARAMETER :: stashcode_pws_cloudturb_pot_diag = 45
INTEGER, PARAMETER :: stashcode_pws_wafc_caturb    =  47
INTEGER, PARAMETER :: stashcode_pws_dustconc_surf  =  58
INTEGER, PARAMETER :: stashcode_pws_dustconc_5000  =  59
INTEGER, PARAMETER :: stashcode_pws_zenithdelay    =  60
INTEGER, PARAMETER :: stashcode_pws_isotherm_ms70_ht    =  61
INTEGER, PARAMETER :: stashcode_pws_isotherm_ms70_press =  62
INTEGER, PARAMETER :: stashcode_pws_isotherm_ms70_icao  =  63
INTEGER, PARAMETER :: stashcode_pws_1p5m_vis_tot        =  70
INTEGER, PARAMETER :: stashcode_pws_1p5m_vis_dust       =  71
INTEGER, PARAMETER :: stashcode_pws_panofsky_turb       =  74
INTEGER, PARAMETER :: stashcode_pws_ellrodt1_turb       =  75
INTEGER, PARAMETER :: stashcode_pws_inv_richardson      =  76
INTEGER, PARAMETER :: stashcode_pws_upd_helicity_5k     =  80
INTEGER, PARAMETER :: stashcode_pws_tropopause_press    =  84
INTEGER, PARAMETER :: stashcode_pws_tropopause_temp     =  85
INTEGER, PARAMETER :: stashcode_pws_tropopause_ht       =  86
INTEGER, PARAMETER :: stashcode_pws_tropopause_icao     =  87
! ---------------------------------------------------------
! Stashcodes which PWS diags depend on
! ---------------------------------------------------------
! Boundary layer, section 3
INTEGER, PARAMETER :: stashcode_bl_windspeed_10mb = 227
INTEGER, PARAMETER :: stashcode_bl_1p5m_temp      = 236
INTEGER, PARAMETER :: stashcode_bl_1p5m_vis_tot   = 281
INTEGER, PARAMETER :: stashcode_bl_tke            = 473
! Convection, section 5
INTEGER, PARAMETER :: stashcode_conv_icao_base    = 210
INTEGER, PARAMETER :: stashcode_conv_icao_top     = 211
! Gravity wave drag, section 6
INTEGER, PARAMETER :: stashcode_gwd_stress_lev_u  = 201
INTEGER, PARAMETER :: stashcode_gwd_stress_lev_v  = 202
! Processed dynamics, section 15
INTEGER, PARAMETER :: stashcode_dyn_wind_ub       = 201
INTEGER, PARAMETER :: stashcode_dyn_wind_vb       = 202
! Processed physics, section 16
INTEGER, PARAMETER :: stashcode_phy_geopht        = 202
! Aerosols
INTEGER, PARAMETER :: stashcode_aero_total_dust   = 257

END MODULE um_stashcode_mod
