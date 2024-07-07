! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    MODULE ACP_NAMEL_MOD  -------------------------------
!
!    Purpose : Read in AC Namelist (&ACP) and process.
!
!   ACP_NAMEL:
!   Set defaults for ACP namelist variables.
!                    Read in namelist and process.
!
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Logical system components covered:
!
!    Project Task : P3
!
!    External documentation:
!
!
!
!    Arguments:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation

!
! Documentation of AC scheme variables (from deleted include files)
! -----------------------------------------------------------------
!   The acp namelist contains parameters controlling the assimilation.
!
!   AC_OBS_TYPES      - List of AC Obs Types to be processed.
!                     - Order of processing is stored in LACT.
!   AC_ORDER          - Order which AC Obs Types MUST be processed.
!                     - (Coded as group*1000+Obstype)
!   CSCFACT_H         - Horizontal Correlation Scale Factor.
!                     - Varies with latitude.
!   CSCFACT_V         - Vertical Correlation Scale Factor.
!                     - Varies with level
!                     - NSLABS_SCFACT(J-1)+1 TO NSLABS_SCFACT(J)
!   DEF_AC_ORDER      - Default Order and groupings of Obs Types.
!   DEF_CSCALE_START  - Default Correlation Scale at start of
!   DEF_CSCALE_OBTIME - insertion period/observation time/end of
!   DEF_CSCALE_END    - insertion period for each group. Use
!                     - CSCALE_START/OBTIME/END to change defaults.
!   DEF_MODE_HANAL    - Default Mode of Horizontal Analysis.
!   DEF_FI_VAR_FACTOR - Default group dep scaling in FI
!   DEF_NO_ANAL_LEVS  - Default No of Analysis Levels for each group.
!   DEF_NO_WT_LEVS    - Default No of Weight Levels for each group.
!                    - Use N_ANAL_LEVS/N_WT_LEVS to change defaults.
!   DEF_NO_ITERATIONS - Default No of Iterations for groups. Use
!                     - NO_ITERATIONS in namelist to change defaults.
!   DEF_INTERVAL_ITER - Default No of Iterations for groups. Use
!                     - INTERVAL_ITER in namelist to change defaults.
!   DEF_OBTHIN        - Default ob thinning (use OBTHIN in namelist)
!                     - (values of 1 imply no thinning, N implies
!                     - 1/N reports assimilated every Nth step)
!   DEF_NUDGE_NH      - Default Nudging Coeffs for NH for groups.
!   DEF_NUDGE_TR      - Default Nudging Coeffs for TR for groups.
!   DEF_NUDGE_SH      - Default Nudging Coeffs for SH for groups.
!                    - Use NUDGE_NH/TR/SH in namelist to change
!                    - defaults.
!   DEF_NUDGE_LAM     - Default Nudging Coeffs for LAM for groups.
!                    - Use NUDGE_LAM in namelist to change defaults.
!   DEF_RADINF        - Default Max Normalised Influence Radius.
!                     - Use RADINF in namelist to change defaults.
!   DEF_TGETOBB )     - Default Time Window before/after obs time to
!   DEF_TGETOBA )     - fetch obs for groups. Use TGETOBB/TGETOBA in
!                     - namelist to change deafults.
!   DEF_TIMEB )       - Default Insertion Period before/after obs time
!   DEF_TIMEA )       - for groups. Use TIMEB/TIMEA in namelist
!                     - to change defaults.
!   DF_COEFF          - Coefficient for DIVFILT
!   DF_SCALE          - DIVFILT scale (metres)
!   DF_SCALE_LEV      - DIVFILT scale for each level
!   DIAG_RDOBS        - Diagnostic Control for subroutine RDOBS.
!   EPSILON_LHN       - Epsilon value for use in LHN
!   F1_506           \                                                 .
!   F2_506            } Parameters for 506 ob weight evaluation
!   F3_506           /
!   ALPHA_LHN         - Alpha value for use in LHN
!   LHN_LIMIT         - Limit on + or - Theta incr rate in LHN (K/day)
!   FI_SCALE_LHN      - Recursive filter scale in m
!   NPASS_RF_LHN      - Number of passes through filter
!   FI_SCALE          - FI (Filtered Increments) scale (metres)
!   FI_SCALE_FACTOR   - FI Scale Factor
!   GEOWT_H           - Latitude weights for Geostrophic Increments.
!   GEOWT_V           - Vertical weights for Geostrophic Increments.
!   GROUP_NO          - Group No of each obs type in LACT.
!   GROUP_FIRST       - Position in LACT of first type of each group.
!   GROUP_LAST        - Position in LACT of last  type of each group.
!   GROUP_INDEX       - Corresponding group in DEF_AC_ORDER for
!                     - groups in GROUP_NO.
!   IOMITOBS          - List of Observations not to be assimilated.
!                     - Use Model Observation Numbers to omit obs.
!   IUNITNO           - Unit No of Cache file to store obs
!                     - between timesteps.
!   L_506_OBERR       - Logical switch to control 506 ob weight calc'n
!   L_LHN             - Logical switch to perform latent heat nudging
!   L_LHN_1A          - Logical switch to perform 1A version of LHN
!   L_LHN_SCALE       - Logical switch to control scaling within LHN
!   L_LHN_SEARCH      - Logical switch to control use of LHN_SEARCH
!   L_VERIF_RANGE     - Logical switch to control verification range
!   L_LHN_LIMIT       - Logical switch to control limiting of increments
!   L_LHN_FACT        - Logical switch to control limiting by 1/alpha
!   L_LHN_FILT        - Logical switch to control filtering of incrs
!   LACT              - List of Obs Types to be processed in order
!                     - of processing.
!   LAC_UARS          - Logical switch for UARS assimilation.
!   LAC_MES           - Logical switch for Mesoscale assimilation.
!   LCHECK_GRID       - Logical switch to control CHECK_OBS call
!   LENACT            - No of obs for each type in group fetched
!                     - this timestep.
!   LGEO  )           - Logical switches to calculate
!   LHN_DIAG          - Logical switch for detailed LHN diagnostics
!   LHN_RANGE         - Max search radius used in LHN_SEARCH. Do not set
!                       to zero, set L_LHN_SEARCH=.FALSE. instead.
!   LHYDR )           - Geostrophic/Hydrstatic Increments.
!   LHYDROL           - Logical switch to calc hydrology incrs.
!   LRADAR            - Logical array to determine which radars to use
!   L_LATLON_PRVER    - Logical switch to verify precip in lat/lon area
!     NORTHLAT        - Co-ords in
!     SOUTHLAT        -            real lat/lon
!     WESTLON         -                         for rain
!     EASTLON         -                                  verif area.
!   L_MOPS_EQUALS_RH  - If .TRUE. then MOPS cloud obs are
!                     - rh values (%), else cloud fractions
!   L_OBS_CHECK       - If .FALSE. then skip check to see if there
!                     - are any obs to assimilate (non-oper run only)
!   LSYN              - Logical switch for Synoptic Insertion.
!   LWBAL_SF          - Controls use of WINDBAL routine for surface wind
!   LWBAL_UA          - Controls use of WINDBAL routine for uair wind
!   MASTER_AC_TYPES   - Master list of AC Observation Types
!                     - known to AC Scheme. See DEF_TYPE.
!   MACDIAG           - Diagnostics control for each timestep.
!                     - See AC on use of this array.
!   MDATADFN          - Mode for Data Density Formula.
!   MGEOWT            - Mode of Latitude Weighting for
!                     - Geostrophic Increments.
!   MGLOSSFN          - GLOSS processing Function Option. (see VANLASS)
!   MHCORFN           - Correlation Function Option.
!   MODEACP           - No of timesteps allowed in dimensioning of
!                     - MACDIAG. Code loops back to MACDIAG(1) if
!                     - TIMESTEP_NO > MODEACP
!   MVINT205          - Options for vertical interp    (LASS/GLOSS  )
!   MRAMPFN           - Mode for Time ramp in HORINF
!   MWTFN             - Mode for Weights Formula.
!   NACT              - No of obs types in LACT.
!   N_GROUPS          - No of groups being processed.
!   NO_OBS_FILES      - No of observation files to be used.
!   NO_SCFACT         - List of obs types on which correlation
!                     - scale factor is not to be applied.
!   NON_DIV_COR       - Factor for Non-Divergent Correction.
!   NON_DIV_COR_10M - As NON_DIV_COR but for 10m wind data
!   NPASS_RF          - Number of passes in RFILT
!   NPROG             - Number set by individual programmers
!                     - doing test work. Numbers in use : 1001 (DR)
!   NRADARS           - No of radars in network
!   NSLABS_SCFACT     - Slab for each level
!                     - (Slab is group of levels with same CSCFACT_V)
!   OBS_FORMAT        - Format of AC Obs file (=2, only one format)
!   OBS_UNITNO        - Unit No of first AC Obs file (=70)
!   OBTIME_NOM        - Nominal Observation Time for Synoptic Insertion
!                     - Mode. Relative time from start of assimilation.
!   RADAR_LAT         - Coordinates of radars
!   RADAR_LON         -     "        "    "
!   RADAR_RANGE       - Max range (km) of reliable radar rain rates
!   RADAR_RANGE_MAX   - Max. range of radar data used for LHN (km)
!   RELAX_CF_LHN      - Relaxation coef used for theta incrs in LHN
!   REMOVE_NEG_LH - Logical switch for removing -ve LH values
!                                       Mark Dixon 23/02/06
!   SPEED_LIMIT305    - Min speed of scatwinds for which observed
!                     - direction is used. (below limit speed only
!                     - is assimilated.
!   THRESH_DL         - threshold (mm/hr) between dry/light rain
!   THRESH_LM         - threshold (mm/hr) between light/moderate rain
!   THRESH_DL         - threshold (mm/hr) between moderate/heavy rain
!   THRESH_RMSF       - threshold (mm/hr) for calcn of rms factor score
!   TIMEF_START       - Start    ) Values for
!   TIMEF_OBTIME      - Obs time ) Time Factor
!   TIMEF_END         - End      ) Ramp Function
!   TROPLAT           - Latitude at which parameters start to change
!                     - from their mid-latitude to tropical values.
!   TROPINT           - Interval over which parameters start to change
!                     - from their mid-latitude to tropical values.
!   TYPE_INDEX        - Position of obs types in LACT in MASTER_AC_TYPES
!   USE_CONV_IN_MOPS  - Logical switch for using convection in MOPS
!   VERT_CUTOFF_SL    - No of scale hts over which single level obs used
!   VERT_CUTOFF_BW    - as VERT_CUTOFF_SL but for bogus winds
!   VERT_CUTOFF_BH    - as VERT_CUTOFF_SL but for bogus humidity
!   VERT_COR_AERO     - Vertical correlation scale for aerosol incrs
!   VERT_COR_SCALE    - Vertical Correlation Scale Coefficient.
!                     - calculated from ACP namelist array CSCALE_VERT
!                       (n,1)=extra tropical temps,(n,2)=tropical temps
!                       (n,3)=extra tropical winds,(n,4)=tropical winds
!   VERT_FILT         - Vertical Filtering of Increments
!                     - from soundings.
!   WB_THETA_UA       - If T WINDBAL will calc theta incs from upper air
!                     - winds.
!   WB_THETA_SF       - If T WINDBAL will calc theta incs from surface
!                     - winds.
!   WB_LAT_CC         - horizontal correlation coeff for WINDBAL
!   WB_VERT_V         - vertical variation in WINDBAL correlations
!   WB_LAND_FACTOR    - extra scaling factor of WINDBAL inc over land
!   WB_LAND_SCALE     - apply WB_LAND_FACTOR scaling if true
!   WB_LonOffset      -) These define a subset of a limited area model
!   WB_LonPts         -) within which WINDBAL will work. They exist to
!   WB_LatOffset      -) select a region on which a multigrid Poisson
!   WB_LatPts         -) solver can be used efficiently. Offsets are
!                     -) from start of full LAM grid (so offsets of
!                     -) zero mean no offset). WB_LonPts & WB_LatPts
!                     -) define the length of the subset in their
!                     -) respective directions.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!   Variables to define the Analysis Grid. Defined in ac_control_mod.
!
!   AG = Analysis Grid
!   MG = Model    Grid
!
!   DEF_AGRES_ROWS - Default Ratio of MG Rows to AG Rows for each group.
!                  - Use AGRES_ROWS in namelist to change defaults.
!   DEF_AGRES_PTS  - Default Ratio of MG Points to AG Points for each
!                  - group. (Equatorwards from AGLATDEC).
!                  - Use AGRES_ROWS and AGRES_PTS in namelist to
!                  - change defaults.
!
!   COLATITUDES & LONGITUDES are in RADIANS, not DEGREES.
!
!   Variables defining Analysis Grid.
!   ---------------------------------
!
!   NROWSAG        - No of rows in AG.
!   NPTSAG         - No of points in each AG row.
!   NPTS0AG        - Offset to first point in each row in AG.
!   NPTSAGMX       - Maximum no of points in AG rows.
!   NPTSAGMN       - Minimum no of points in AG rows.
!   MIN_AGPTS      - Minimum no of points on any AG ROW.
!   LAGNP          - Logical set to T if first row of AG = North Pole
!   LAGSP          - Logical set to T if last  row of AG = South Pole
!   STAGROW1       - Stagger of first row of MG to first row of AG,
!                  - (expressed as fraction of AG row spacing).
!   STAGPT1        - Stagger of first point of MG to first point of AG,
!                  - (expressed as fraction of AG row spacing).
!   ROW1AG         - Co-latitude of first row   of AG.
!   PT1AG          - Longitude   of first point of AG.
!   DLATAG         - Co-latitude of row spacing of AG.
!   DLONGAG        - Longitude spacing for each AG row.
!   AGLATDEC       - Co-latitude at which no of pts in AG rows start
!                  - to decrease.
!   AGROWLEN       - Length of row in radians of latitude.
!   COSROWAG       - Array storing 1/(DLONGAG(x)*cos(lat of row x))
!                  - where x is the AG row no.
!
!   Variables defining Model Grid.
!   ------------------------------
!
!   ROW1MG         - Co-latitude of first row of MG in use.
!   ROW1MGTH       - Co-latitude of first row of MG for Theta.
!   ROW1MGUV       - Co-latitude of first row of MG for Winds.
!   PT1MGTH        - Longitude of first point on MG for Theta.
!   PT1MGUV        - Longitude of first point on MG for Winds.
!   DLATMG         - Co-latitude of row spacing of MG (+ve for N->S)
!   DLONGMG        - Longitude spacing of MG.
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!     CONTAINS PARAMETERS CONTROLLING DIAGNOSTIC PRINTOUT.
!     DEFINED IN AC_CONTROL_MOD.
!
!     LDIAGAC =.FALSE. SWITCHES OFF ALL DIAGNOSTICS
!     LLDAC(#)=.FALSE. SWITCHES OFF DIAGNOSTICS OF TYPE #
!
!     DIAGNOSTIC TYPES & THEIR ASSOCIATED PARAMETERS ARE :-
!     --------------------------------------------------
!
!     LLDAC(1) = .TRUE.
!     -----------------
!
!     DETAILED DIAGNOSTICS ON A SHORT LIST OF OBSERVATIONS
!     TO BE DONE FROM AC ETC.
!     THIS TYPE OF DIAGNOSTIC IS SET-UP & USED BY TYPE 5 IF LLDAG0=T.
!     AFTER EACH ITERATION THE LIST OF OBS INFLUENCING THE SELECTED PNT
!     IS ADDED TO MDACO. THUS FOR SUBSEQUENT ITERATIONS DETAILED
!     DIAGNOSTICS ON THESE OBS ARE OBTAINED. TO GET THIS MODE SET
!     LLDAC(5)=T LLDAG0=T.
!
!     THIS TYPE CAN ALSO BE SWITCHED ON INDEPENDENTLY OF TYPE 5 BY :-
!     MODACO  MODE FOR SETTING UP LIST MDAC:-
!             1 TAKE FIRST NDACP OBS IN CURRENT VECTOR
!             2 SEARCH FOR THOSE OBS IN MDACO DEFINED IN NAMELIST ADIAG
!             3 SEARCH FOR THOSE OBS IN LTD AREA FROM    NAMELIST ADIAG
!
!     NDACOP  MAX NO OF OBS ON WHICH DIAGNOSTICS ARE REQUIRED
!     NDACO       NO OF OBS ON WHICH DIAGNOSTICS ARE REQUIRED
!     MDACO     LIST OF OBS ON WHICH DIAGNOSTICS ARE REQUIRED
!               POINTS TO POSITION OF OB IN COMOBS.
!
!     NDACP  MAX NO OF OBS OF CURRENT TYPE ON WHICH DIAGS REQD.
!     NDAC       NO OF OBS OF CURRENT TYPE ON WHICH DIAGS REQD.
!     MDAC     LIST OF OBS OF CURRENT TYPE ON WHICH DIAGS REQD.
!               POINTS TO POSITION OF OB IN CURRENT VECTORS FROM GETOBS.
!
!     NDACVP  MAX NO OF PARAMETERS WHICH CAN BE STORED FOR EACH OB
!     NDACV       NO OF PARAMETERS WHICH ARE    STORED FOR EACH OB
!     DACV        STORED PARAMETERS FOR EACH OB
!     CACV        DESCRIPTION OF EACH PARAMETER
!     CACT1/2/3/4 TITLES SET UP BY DACOIN
!     LTITLE      CONTROLS TITLING OF OUTPUT FROM DACOUT
!
!     LLDAC(2) = .TRUE.
!     -----------------
!     STATISTICS OF OBSERVATION-MODEL INCREMENTS
!     PRINTED OUT IN DIAGO
!     MDIAGTPS: LIST OF TYPES TO BE PROCESSED BY DIAGO.
!               SET MDIAGTPS(1)=0 FOR 'ALL'(DEFAULT).
!     MDIAG:    1 = ONLY CALLS FROM AC PROCESSED
!               2 = ONLY CALLS FROM Van### BEFORE VRTF PROCESSED
!               4 = ONLY CALLS FROM Van### AFTER VRTF PROCESSED
!                   (Van### is group of vertical analysis routines)
!               0 = ALL CALLS PROCESSED (DEFAULT)
!               BINARY COMBINATIONS SUPPORTED.
!     LLBAND:   F = GLOBAL STATISTICS (DEFAULT).
!               T = SEPARATE STATISTICS FOR BANDS 90N-22N,
!                                                 22N-22S,
!                                                 22S-90S.
!     LRMS   :   T/F = Print RMS/Mean Square Values. Default = T.
!     LNORMF :   T : Use Normalisation Factors (NF) as weights (Default)
!            :   F : Set NF to 1.0 (ie. no weights used).
!     LTEMP  :   T : Print Temperature statistics. Default. Theta
!                    Increments are converted to temperature increments
!                    and p* is assumed to be 1000mb for all obs.
!            :   F : Print Theta statistics.
!     LVERIF :   T : Sets parameters to get statistics for verification
!                    purposes. LRMS=T,LTEMP=T,LNORMF=F,LGEO=F,LHYDR=F
!            :   F : Default.
!
!     LLDAC(3) = .TRUE.
!     -----------------
!     STATISTICS OF INCREMENTS ON ANALYSIS GRID
!
!     LLDAC(4) = .TRUE.
!     -----------------
!     STATISTICS OF INCREMENTS ON MODEL GRID
!     PRINTED OUT IN MMSPT
!     MEAN AND MEAN-SQUARE OF INCREMENTS PRINTED FOR :-
!         1. THE WHOLE MODEL(NOT WITH GEOSTROPHIC INCREMENTS)
!     AND 2. EACH LEVEL OF THE MODEL
!
!     LLDAC(5) = .TRUE.
!     -----------------
!     DETAILS OF OBS INFLUENCING A SPECIFIED POINT.
!     PARAMETERS IN NAMELIST ADIAG :-
!     DAGLAT  -  LATITUDE  OF SPECIFIED POINT. DEFAULT =  57.0
!     DAGLON  -  LONGITUDE OF SPECIFIED POINT. DEFAULT = 340.0
!     LLDAG0      = SWITCH ON OPTION 1 FOR THOSE OBS FOUND, SO THAT
!                    THEIR DETAILS ARE PRINTED NEXT TIME THROUGH AC.
!     LLDAG(NAVT) = DO DIAGNOSTIC  FOR THIS VARIABLE TYPE.
!
!-----------------------------------------------------------------------


MODULE acp_namel_mod

USE comobs_mod, ONLY: nobtypmx
USE ac_control_mod
USE missing_data_mod, ONLY: imdi, rmdi
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Switches for assm mode
LOGICAL :: l_ac = .TRUE.

! Time at which data assimilation starts (Hours after Basis Time)
INTEGER :: a_assim_start_min = imdi

! Time at which data assimilation  ends
INTEGER :: a_assim_end_min   = imdi

INTEGER :: ac_order (nobtypmx)
INTEGER :: timeb    (nobtypmx)
INTEGER :: timea    (nobtypmx)
INTEGER :: tgetobb  (nobtypmx)
INTEGER :: tgetoba  (nobtypmx)
REAL    :: nudge_lam(nobtypmx)

NAMELIST /acp/                                                   &
ac_order, ac_obs_types, lac_mes, no_obs_files, obs_format,       &
l_lhn, l_lhn_1a, lhn_range, l_lhn_scale, l_lhn_search, lhn_diag, &
epsilon_lhn, alpha_lhn, relax_cf_lhn, l_lhn_limit, lhn_limit,    &
l_lhn_fact, l_lhn_filt, fi_scale_lhn, npass_rf_lhn,              &
remove_neg_lh, nudge_lam, diag_rdobs, macdiag, use_conv_in_mops, &
timea, timeb, tgetoba, tgetobb, a_assim_start_min, a_assim_end_min

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACP_NAMEL_MOD'

CONTAINS

SUBROUTINE acp_namel (bl_levels, tr_levels,            &
                      p_rows, u_rows,row_len,          &
                      timestep, icode, cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore,    ONLY: mype
USE Field_Types,   ONLY: fld_type_p
USE def_group_mod, ONLY: def_group
USE def_type_mod, ONLY: def_type
USE group_dep_var_mod, ONLY: group_dep_var
USE type_dep_var_mod, ONLY: type_dep_var

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE check_iostat_mod
USE ereport_mod, ONLY: ereport

USE model_domain_mod, ONLY: model_type, mt_global
USE nlsizes_namelist_mod, ONLY: model_levels

IMPLICIT NONE


! Imported variables (INTENT=IN):L
INTEGER ::  bl_levels           ! total number of boundary layer
                                ! levels
INTEGER ::  tr_levels           ! total number of tracer levels
INTEGER ::  p_rows              ! Number of rows (for pstar)
INTEGER ::  u_rows              ! Number of rows (for wind)
INTEGER ::  row_len             ! Length of each row.

REAL ::     timestep            ! timestep in seconds

! Exported variables (INTENT=OUT)
INTEGER :: icode ! Non zero for failure

CHARACTER(LEN=errormessagelength) :: cmessage          ! Reason for failure
!
! ---------------------------------------------------------------------
INTEGER :: p_rows_global,u_rows_global
! ---------------------------------------------------------------------
REAL :: geowt_nh
REAL :: geowt_sh
REAL :: nudge_nh(nobtypmx)
REAL :: nudge_tr(nobtypmx)
REAL :: nudge_sh(nobtypmx)
REAL :: scfact_nh
REAL :: scfact_tr
REAL :: scfact_sh
REAL :: wb_cc_n                     ! Northern hemisphere
                                    ! corln. coeff. for WINDBAL.
REAL :: wb_cc_s                     ! Southern hemisphere
                                    ! corln. coeff. for WINDBAL.
REAL :: geowt_lam
REAL :: wb_cc_lam                   ! WINDBAL correlation

REAL :: geowt_vert(21)
REAL :: cscale_vert(21,4)
REAL :: cscale_vert_g(21,4)
REAL :: def_cscale_vert_mes(21,4)
REAL :: scfact_vert(model_levels_max)
REAL :: fi_var_factor(nobtypmx)
REAL :: temp_coord
REAL :: cscale_vert_aero

INTEGER :: no_iterations(nobtypmx)
INTEGER :: interval_iter(nobtypmx)
INTEGER :: n_anal_levs  (nobtypmx)
INTEGER :: n_wt_levs    (nobtypmx)
INTEGER :: agres_rows   (nobtypmx)
INTEGER :: agres_pts    (nobtypmx)
INTEGER :: mode_hanal   (nobtypmx)
INTEGER :: radinf       (nobtypmx)
INTEGER :: obthin       (nobtypmx)
INTEGER :: cscale_start (nobtypmx)
INTEGER :: cscale_obtime(nobtypmx)
INTEGER :: cscale_end   (nobtypmx)
INTEGER :: j,jj,jlev,jrow,jtyp,jobt
INTEGER :: ncount
INTEGER :: wb_start_lev                 ! First level for non zero
                                     ! vertical correlation coeff.
INTEGER :: wb_end_lev                   ! Last  level for non zero
                                     ! vertical correlation coeff.

INTEGER :: errorstatus
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'ACP_NAMEL'
CHARACTER (LEN=errormessagelength) :: iomessage
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! Vertical correlation scale defaults; 21 values at 50hPa intervals
! for extratropics T,tropics T,for extratropics V,tropics V
DATA cscale_vert_g/                                                 &
       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.5, 5.0, 4.5, 4.0,         &
        3.5, 3.0, 3.0, 3.0, 3.0, 4.0, 3.5, 3.0, 2.5, 1.5, 1.5,    &
        6.0, 6.0, 5.0, 5.0, 5.0, 5.0, 4.0, 3.5, 3.5, 3.0,         &
        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5,    &
        6.5, 6.5, 6.0, 6.0, 5.5, 5.0, 5.0, 4.0, 3.5, 3.5,         &
        3.5, 3.5, 3.0, 3.0, 3.0, 3.0, 3.0, 2.5, 2.0, 1.5, 1.5/
! for LAM Tropics and Extra-tropics set the same
DATA cscale_vert/                                                 &
       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
        9.0, 6.0, 5.0, 5.0, 5.0, 5.0, 4.0, 3.5, 3.5, 3.0,         &
        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5,    &
        9.0, 6.0, 5.0, 5.0, 5.0, 5.0, 4.0, 3.5, 3.5, 3.0,         &
        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5/

DATA def_cscale_vert_mes/                                         &
       18.0,18.0,18.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
       18.0,18.0,18.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
       27.0,27.0,27.0,12.0, 8.0, 6.0, 4.0, 3.5, 3.5, 3.0,         &
        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5,    &
       27.0,27.0,27.0,12.0, 8.0, 6.0, 4.0, 3.5, 3.5, 3.0,         &
        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5/

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

p_rows_global=glsize(2,fld_type_p)
u_rows_global=glsize(2,fld_type_p)-1

! Hard-wired settings
! Switch  for UARS assimilation.
lac_uars = .FALSE.

! Switch for 'non-oper no obs' test in AC1
l_obs_check = .TRUE.

DO jobt = 1, nobtypmx
  ! Group dependent arrays

  no_iterations(jobt)  = imdi
  interval_iter(jobt)  = imdi
  n_anal_levs(jobt)    = imdi
  n_wt_levs(jobt)      = imdi
  agres_rows(jobt)     = imdi
  agres_pts (jobt)     = imdi
  mode_hanal(jobt)     = imdi
  fi_var_factor(jobt)  = rmdi

  IF (model_type == mt_global) THEN
    nudge_nh(jobt)       = rmdi
    nudge_tr(jobt)       = rmdi
    nudge_sh(jobt)       = rmdi
  END IF

  ! Obs Type dependent arrays in namelist.
  radinf(jobt)         = 0
  obthin(jobt)         = 0
  cscale_start(jobt)   = 0
  cscale_obtime(jobt)  = 0
  cscale_end(jobt)     = 0
  no_scfact(jobt)      = 0

END DO

no_scfact(1) = 305

IF (model_type == mt_global) THEN

  ! Global version:
  ! Nominal ob time relative to start of assm for synoptic insertion
  obtime_nom = 300.0

  ! Correlation Scale Factors (Horizontal)
  IF (lac_uars) THEN
    scfact_nh = 1.0
    scfact_tr = 1.0
    scfact_sh = 1.0

  ELSE
    scfact_nh = 1.0
    scfact_tr = 1.3
    scfact_sh = 1.5

  END IF

  ! Scaling for geostrophic increments
  IF (lac_uars) THEN
    geowt_nh  = 0.7
    geowt_sh  = 0.7

  ELSE
    geowt_nh  = 0.5
    geowt_sh  = 0.7

  END IF

  ! Latitude at which mid-latitude Analysis Parameters
  ! begin to change to tropical values.
  troplat=30.0

  ! Interval (in degrees Lat) over which Relaxation Coefficients
  ! and Correlation Scale Factor change.
  tropint = 10.0

  ! Latitude at which model and Analysis Grid become identical.
  ! analysis grid=model grid in MPP mode
  aglatdec = 90.0

ELSE
  ! Limited Area Version

  ! Nominal ob time relative to start of assm for synoptic insertion
  obtime_nom = 150.0

  ! Scaling for geostrophic increments
  geowt_lam = 0.5

  ! Force Model and Analysis Grid to be identical for Limited Area
  aglatdec = 90.0
END IF

! Parameters for time factor
IF (lac_uars) THEN
  timef_start  = 0.01
  timef_end    = 0.01

ELSE
  timef_start  = 0.05
  timef_end    = 0.05

END IF

timef_obtime = 1.0

! Min no. of points on analysis grid row
min_agpts = 14

! Geostrophic & hydrostatic incs.
lgeo  = .TRUE.
lhydr = .TRUE.

! Hydrology incs.
lhydrol = .TRUE.

! Variable (rh or cloud fraction) used in MOPS obs
l_mops_equals_rh = .FALSE.

! Lat/lon area selection for precip verif
l_latlon_prver = .FALSE.
! Default lat/lon area when above logical set to .TRUE.
! (this area is larger than current FRONTIERS area)
northlat = 60.0
southlat = 45.0
westlon  = -20.0
eastlon  = 10.0

! Balanced surface pressure & pot. temp. increments for wind ?
IF (model_type == mt_global) THEN
  lwbal_sf = .TRUE.
  lwbal_ua = .TRUE.

ELSE
  lwbal_sf = .FALSE.
  lwbal_ua = .FALSE.

END IF

! Enable calculation of balanced pot. temp. incs. for wind ?
! if true WINDBAL theta incrs are calculated from upper air winds.
wb_theta_ua = .TRUE.

! if true WINDBAL theta incrs are calculated from surface winds
! (using HYDRST).
wb_theta_sf = .TRUE.

! Initialize variables to define the correlation coefficients for
! WINDBAL.
IF (model_type == mt_global) THEN
  wb_cc_n       = 0.70
  wb_cc_s       = 0.80
ELSE
  wb_cc_lam     = 0.70

  ! Initialize variables to define a subset of limited area model on
  ! which a multigrid can be run efficiently (for use in WINDBAL).
  ! The LAM values are given (with mesoscale values after the comment).
  ! N.B. the LAM is 229 E/W by 132 N/S, and the mesoscale model is 92 by
  ! 92.

END IF
wb_start_lev   = bl_levels - 1     !first non zero level
wb_end_lev     = model_levels - 2  !last  non zero level
wb_land_scale  = .TRUE.  ! rescale WINDBAL incs over land ?
wb_land_factor = 0.5     ! factor for rescaling   "   "

! Synoptic insertion mode
lsyn = .FALSE.

! Vertical filtering parameter (used in vrtf)
vert_filt = 0.0

! vertical correlation function coeff for aerosol
cscale_vert_aero = 18.0

! Vertical correlation function coeff
! EXP(-(CSCALE_VERT*(LN(P1)-LN(P2)))**2)
! VERT_COR_SCALE is interpolated from  CSCALE_VERT
! which is defined for a regular set of 21 eta levels
! (1 through 0 step 0.05)
! 4 copies for extratropical non-wind (,1);tropical non-wind (,2)
!          for extratropical     wind (,3);tropical     wind (,4)
IF (lac_uars) THEN
  ! Altered to decay more sharply than before
  ! (b in equation 5.1.3 is 9 rather than 3 (old default)
  ! or 16 in tropics - same at all levels)
  DO j = 1, 21
    cscale_vert(j,1) = 3.0
    cscale_vert(j,2) = 4.0
    cscale_vert(j,3) = 3.0
    cscale_vert(j,4) = 4.0
    cscale_vert_g(j,1) = 3.0
    cscale_vert_g(j,2) = 4.0
    cscale_vert_g(j,3) = 3.0
    cscale_vert_g(j,4) = 4.0
  END DO
ELSE IF (lac_mes) THEN
  ! use revised mes values
  DO jj = 1, 4
    DO j = 1, 21
      cscale_vert(j,jj) = def_cscale_vert_mes(j,jj)
      cscale_vert_g(j,jj) = def_cscale_vert_mes(j,jj)

    END DO
  END DO
  !     ELSE for LAM and GL (oper) use DATA table
END IF

! Cut off in vertical correlation function
! single lvl obs(203,303,403) infl only vert_cutoff_sl scale hts
! same cut-off correl. applies to extrap'n of incomplete soundings.
IF (lac_uars) THEN
  vert_cutoff_sl = 1.0

ELSE
  vert_cutoff_sl = 2.0

END IF

! Cut off in vertical correlation function for bogus wind data
IF (lac_uars) THEN
  vert_cutoff_bw = 1.0

ELSE
  vert_cutoff_bw = 1.0

END IF

! Cut off in vertical correlation function for bogus humidity data
vert_cutoff_bh = 0.10

! Correlation function option
!     MHCORFN = 3 For auto-regressive function
!             = 4 with mod which tends to zero at edge on infl
!               +16 For time factor option in INITAC & AC
IF (lac_mes) THEN
  mhcorfn = 3

ELSE
  mhcorfn = 4

END IF

! Set ac diagnostic mode (each timestep)
! +8 FOR DIAG ONLY ITERATION AT END OF ITERATIONS
! +32 FOR DIAGS ON (OVERRULING LDIAGAC)
! +64 FOR DIAG ONLY ITERATION BETWEEN ITERATIONS

IF (lac_mes) THEN
  ! Set MACDIAG to get diagnostics on timesteps 1, 12, 24, 36.
  macdiag(1)   = 32
  macdiag(12)  = 32
  macdiag(24)  = 32
  macdiag(36)  = 32

ELSE  !  global and lam
  ! Set MACDIAG to get diagnostics on timesteps 1, 18 and 36.
  macdiag(1)  = 32
  macdiag(18) = 32
  macdiag(36) = 32

END IF

! set mode for weights formula
mwtfn = 2

! set mode for time ramp (1=linear,2=quadratic)
IF (lac_uars) THEN
  mrampfn = 1

ELSE
  mrampfn = 2

END IF

! set mode for data density formula
mdatadfn = 2

! set mode for GLOSS (type 208) processing option
! 1 as 205,2 use GLOSS error ratio,
! 3 do constrained retrieval, 4 as 3 without a background correction
mglossfn = 1

!  Min speed limit for ERS-1 winds above which direction is believed
speed_limit305 = 4.0

! set parameters to define LASS/GLOSS vertical interpolation
mvint205 = 2     ! 1 for V_INT_T, 2 for V_INT_TP

! set mode for latiude weighting of geostrophic increments
mgeowt = 3

!     VERTICAL SCALING FACTOR FOR GEOSTROPHIC INCREMENTS
!     this is defined for a regular set of 21 eta levels
!       (1 through 0 step 0.05)
IF (lac_uars) THEN
  ! UARS defaults give full weights, except near jet level
  DO j=1,21
    geowt_vert(j) = 1.0
  END DO

  geowt_vert(13)=0.8
  geowt_vert(14)=0.6
  geowt_vert(15)=0.4
  geowt_vert(16)=0.3
  geowt_vert(17)=0.4
  geowt_vert(18)=0.6
  geowt_vert(19)=0.8
ELSE
  ! the defaults here give full wt below 500mb and zero wt above 250mb
  DO j=1,10
    geowt_vert(j) = 1.0
  END DO

  geowt_vert(11) = 0.820
  geowt_vert(12) = 0.720
  geowt_vert(13) = 0.620
  geowt_vert(14) = 0.500
  geowt_vert(15) = 0.367
  geowt_vert(16) = 0.200

  DO j = 17, 21
    geowt_vert(j) = 0.0
  END DO
END IF

! Vertical scaling for horizontal correlation scale
IF (lac_mes) THEN
  DO j = 1, bl_levels
    scfact_vert(j)   = 0.875
    nslabs_scfact(j) = 0
    cscfact_v(j)     = 0.0
  END DO
  DO j = bl_levels+1, model_levels_max
    scfact_vert(j)   = 1.0
    nslabs_scfact(j) = 0
    cscfact_v(j)     = 0.0
  END DO

ELSE
  DO j = 1, model_levels_max
    scfact_vert(j)   = 1.0
    nslabs_scfact(j) = 0
    cscfact_v(j)     = 0.0
  END DO

END IF

! Defaults for Filtered Increment Mode of Horizontal Analysis
! For testing, FI_SCALE_FACTOR=1.0 will match current HORINF
! Eventually, higher levels can be given a larger scale.
! DF_COEFF=1.0 is initial estimate. DF_COEFF=0.0 turns off DIVFILT
! NPASS_RF=2 matches SOAR and is default.

npass_rf = 2
IF (lac_mes) THEN
  !       For MOPS data, use FI method without filtering on model grid.
  !       (Assumes that no other mes data groups require standard FI)
  npass_rf = 0
END IF

! Only MOPS data are assimilated now so that LAC_MES is T and
! NPASS_RF is 0. If this is not the case abort
IF (npass_rf /= 0) THEN

  WRITE(umMessage,*) 'NPASS_RF NOT ZERO'
  CALL umPrint(umMessage,src='acp_namel')
  WRITE(umMessage,'(A)')'Check that MOPS data assimilation is on in gui'
  CALL umPrint(umMessage,src='acp_namel')
  cmessage='MOPS data not switched on in input namelists'
  icode=1
  GO TO 9999
END IF

fi_scale = 400000.0
df_scale = 400000.0

DO jlev = 1, model_levels_max
  fi_scale_factor(jlev) = 1.0
  df_coeff(jlev)        = 1.0
  df_scale_lev(jlev)    = df_scale * fi_scale_factor(jlev)

END DO

! Factor for non-divergence correction term in correlations
non_div_cor = 0.8
!  same factor for 10m wind data
IF (lac_mes) THEN
  non_div_cor_10m = 0.0
ELSE
  non_div_cor_10m = 0.8
END IF

! Thresholds (mm/hr) for rainfall verification
thresh_dl   = 0.03
thresh_lm   = 0.125
thresh_mh   = 0.5
thresh_rmsf = 0.125

! Radar data reliability ranges.
! RADAR_RANGE is limit of high reliability (used in LHN & HCS)
! RADAR_RANGE_MAX is limit of usability (used in LHN)
radar_range = 100.0
!  Max Radar range
radar_range_max = 200.0
!  --- SPECIFY RADAR LOCATIONS ------
radar_lat(1)   =  58.21
radar_lon(1)   =  353.81      !  BEACON HILL

radar_lat(2)   =  57.43
radar_lon(2)   =  357.97      !  HILL OF DUDWICK

radar_lat(3)   =  55.69
radar_lon(3)   =  355.77      !  CORSE HILL

radar_lat(4)   =  54.50
radar_lon(4)   =  353.65      !  CASTOR BAY

radar_lat(5)   =  53.43
radar_lon(5)   =  353.75      !  DUBLIN

radar_lat(6)   =  52.69
radar_lon(6)   =  351.07      !  SHANNON

radar_lat(7)   =  53.75
radar_lon(7)   =  357.72      !  HAMELDON

radar_lat(8)   =  53.33
radar_lon(8)   =  359.45      !  INGHAM

radar_lat(9)   =  52.40
radar_lon(9)   =  357.40      !  CLEE HILL

radar_lat(10)  =  51.69
radar_lon(10)  =  359.47      !  CHENIES

radar_lat(11)  =  51.98
radar_lon(11)  =  355.55      !  CRUGYGORLLWYN

radar_lat(12)  =  50.96
radar_lon(12)  =  356.55      !  COBBACOMBE

radar_lat(13)  =  50.82
radar_lon(13)  =  357.44      !  WARDON HILL

radar_lat(14)  =  50.00
radar_lon(14)  =  354.77      !  PREDANNACK

radar_lat(15)  =  49.18
radar_lon(15)  =  357.78      !  JERSEY

! List of which radars to use
! (.TRUE. implies MOPS precip data close to that radar will be used
! in rainfall verification and hydrology correction)
lradar(1)   = .TRUE.     !  Beacon Hill
lradar(2)   = .TRUE.     !  Hill of Dudwick
lradar(3)   = .TRUE.     !  Corse Hill
lradar(4)   = .TRUE.     !  Castor Bay
lradar(5)   = .TRUE.     !  Dublin
lradar(6)   = .TRUE.     !  Shannon
lradar(7)   = .TRUE.     !  Hameldon
lradar(8)   = .TRUE.     !  Ingham
lradar(9)   = .TRUE.     !  Clee
lradar(10)  = .TRUE.     !  Chenies
lradar(11)  = .TRUE.     !  Crugygorllwyn
lradar(12)  = .TRUE.     !  Cobbacombe
lradar(13)  = .TRUE.     !  Wardon Hill
lradar(14)  = .TRUE.     !  Predannack
lradar(15)  = .TRUE.     !  Jersey


!  506 observation error calculations (WP171, eqn10)
l_506_oberr = .FALSE.
f1_506 = 10.0
f2_506 = 1.0
f3_506 = 0.0

!  Verify up to reliable range, or max range
l_verif_range = .TRUE.

! List of observations to be omitted from assimilation
! ** Specify model observation type to omit obs **
DO jtyp = 1, nanaltyp
  iomitobs(jtyp) = 0

END DO

! Diagnostics options for programmers
nprog = 0

! Set timing switch
ltimer_ac = .FALSE.

! Set control for calling CHECK_OBS
lcheck_grid = .FALSE.

!  2   Check Validity of Namelist Parameters
! Check AGLATDEC
IF (aglatdec <  89.0) THEN
  aglatdec = 90.0
  CALL umPrint(' AGLATDEC reset to 90. when running MPP ', &
      src='acp_namel',pe=0)
END IF
! Check that AC_OBS_TYPES has any obs types
ncount = 0

DO jobt = 1, nobtypmx
  IF (ac_obs_types(jobt)  >   0) THEN
    ncount = ncount + 1

  END IF
END DO

IF (ncount  ==  0) THEN
  icode    = 1
  cmessage = ' ACP_NAMEL : No obs types in AC_OBS_TYPES.'
  GO TO 9999

END IF

! Convert any old type numbers to new type numbers.
DO jobt = 1, nobtypmx
  IF (ac_obs_types(jobt)  ==  501) THEN
    ac_obs_types(jobt) = 302
    CALL umPrint('Type 501 in AC_OBS_TYPES changed to Type 302', &
        src='acp_namel',pe=0)
  END IF

  IF (ac_obs_types(jobt) == 502) THEN
    ac_obs_types(jobt) = 305
    CALL umPrint('Type 502 in AC_OBS_TYPES changed to Type 305', &
        src='acp_namel',pe=0)
  END IF
END DO

IF (model_type == mt_global) THEN
  ! Check that troplat and tropint are correctly specified
  IF (troplat + tropint  >   90.0) THEN
    icode    = 1
    cmessage =                                                      &
            ' ACP_NL : TROPLAT/TROPINT bug, is TROPLAT-TROPINT >= 0'
    GO TO 9999

  END IF

  ! Check geowt_nh non zero if geostrophic increments are to be used
  IF (geowt_nh  ==  0.0 .AND. (lgeo)) THEN
    icode    = 1
    cmessage = 'ACP_NL : GEOWT_NH = 0 not allowed with LGEO = T'
    GO TO 9999

  END IF

ELSE
  ! Check geowt_lam non zero if geostrophic increments are to be used
  IF (geowt_lam  ==  0.0 .AND. (lgeo)) THEN
    icode    = 1
    cmessage = 'ACP_NL : GEOWT_LAM=0 not allowed with LGEO = T'
    GO TO 9999

  END IF
END IF  ! if GLOBAL

! Check windbal variables
IF (model_type == mt_global) THEN
  IF (wb_cc_n  <   0.0 .OR. wb_cc_n  >   1.0) THEN
    icode    = 1
    cmessage = 'ACP_NL : invalid WB_CC_N'
    GO TO 9999

  END IF
  IF (wb_cc_s  <   0.0 .OR. wb_cc_s  >   1.0) THEN
    icode    = 1
    cmessage = 'ACP_NL : invalid WB_CC_S'
    GO TO 9999

  END IF
ELSE
  IF (wb_cc_lam  <   0.0 .OR. wb_cc_lam  >   1.0) THEN
    icode    = 1
    cmessage = 'ACP_NL : invalid WB_CC_LAM'
    GO TO 9999

  END IF
END IF  ! if GLOBAL

! Check MRAMPFN
IF (mrampfn  /=  1 .AND. mrampfn  /=  2) THEN
  icode    = 1
  cmessage = 'ACP_NL : invalid MRAMPFN'
  GO TO 9999

END IF

! Check MDATADFN
IF (mdatadfn  /=  1 .AND. mdatadfn  /=  2) THEN
  icode    = 1
  cmessage = 'ACP_NL : invalid MDATADFN'
  GO TO 9999

END IF

! Check MGLOSSFN
IF (mglossfn  <   1 .OR. mglossfn  >   4) THEN
  icode    = 1
  cmessage = 'ACP_NL : invalid MGLOSSFN'
  GO TO 9999

END IF

! Check MWTFN
IF (mwtfn /= 1 .AND. mwtfn /= 2) THEN
  icode    = 1
  cmessage = 'ACP_NL : invalid MWTFN'
  GO TO 9999

END IF

! Check MHCORFN
IF (MOD(mhcorfn,8) <  1 .OR. MOD(mhcorfn,8) >  4) THEN
  icode    = 1
  cmessage = 'ACP_NL : invalid MHCORFN'
  GO TO 9999

END IF

! Check MVINT205
IF (mvint205  /=  1 .AND. mvint205  /=  2) THEN
  icode    = 1
  cmessage = 'ACP_NL : invalid MVINT205'
  GO TO 9999

END IF

! Check NO_OBS_FILES
IF (no_obs_files >  10) THEN
  icode    = 1
  cmessage = 'ACP_NL : invalid NO_OBS_FILES (max allowed is 10)'
  GO TO 9999
END IF

! Check OBS_FORMAT - only allow portable UM format
IF (obs_format /= 2 .AND. obs_format /= 3) THEN
  icode    = 1
  cmessage = 'ACP_NL : OBS_FORMAT not equal to 2 or 3'
  GO TO 9999
END IF

!  3. Set Defaults for Group Dependent arrays
CALL def_group (bl_levels, tr_levels,         &
 icode, cmessage)

IF (icode >  0) GO TO 9999

!  4. Set Defaults for Obs Type Dependent arrays
CALL def_type (icode,cmessage)

IF (icode >  0) GO TO 9999

!  5. Check/Process Group dependent variables
CALL group_dep_var (ac_order, no_iterations, interval_iter,       &
     n_anal_levs, n_wt_levs,                                      &
     nudge_nh, nudge_tr, nudge_sh,                                &
     nudge_lam,                                                   &
     agres_rows, agres_pts, mode_hanal, fi_var_factor,            &
     icode,cmessage)

IF (icode >  0) GO TO 9999

!  6. Check/Process Obs Type dependent variables
CALL type_dep_var (timeb,timea,tgetobb,tgetoba,radinf,obthin,     &
                   cscale_start,cscale_obtime,cscale_end,         &
                   icode,cmessage)

IF (icode >  0) GO TO 9999


!  11. Convert rainfall thresholds from mm/hr to kgm-2s-1
thresh_dl   = thresh_dl   / 3600
thresh_lm   = thresh_lm   / 3600
thresh_mh   = thresh_mh   / 3600
thresh_rmsf = thresh_rmsf / 3600

!   12. Check Lat/Lon verification coordinates

! Need longitudes in range :  -180 < LON <= 180
IF (eastlon  >   180.0) THEN
  eastlon = eastlon - 360.0
  WRITE(umMessage,*)                                         &
      "EASTLON changed from ",eastlon+360.0," to ",eastlon
  CALL umPrint(umMessage,src='acp_namel',pe=0)
END IF

IF (westlon  >   180.0) THEN
  westlon = westlon - 360.0
  WRITE(umMessage,*)                                         &
      "WESTLON changed from ",westlon+360.0," to ",westlon
  CALL umPrint(umMessage,src='acp_namel',pe=0)
END IF

IF (northlat  <   southlat) THEN
  WRITE(umMessage,*)                                         &
      "Swapping NORTHLAT and SOUTHLAT so that North>South"
  temp_coord = northlat
  northlat   = southlat
  southlat   = temp_coord
  CALL umPrint(umMessage,src='acp_namel',pe=0)
END IF

IF (eastlon  <   westlon) THEN
  WRITE(umMessage,*)                                         &
      "Swapping EASTLON and WESTLON so that East > West"
  temp_coord = eastlon
  eastlon    = westlon
  westlon    = temp_coord
  CALL umPrint(umMessage,src='acp_namel',pe=0)
END IF

9999  CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE acp_namel


SUBROUTINE print_nlist_acp()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer

CALL umPrint('Contents of namelist acp', &
    src='acp_namel')

WRITE(lineBuffer,*)' AC_ORDER = ',ac_order
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' AC_OBS_TYPES = ',ac_obs_types
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' TIMEA = ',timea
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' TIMEB = ',timeb
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' TGETOBA = ',tgetoba
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' TGETOBB = ',tgetobb
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' NPASS_RF = ',npass_rf
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' LAC_MES = ',lac_mes
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' OBS_FORMAT = ',obs_format
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' NO_OBS_FILES = ',no_obs_files
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' DIAG_RDOBS = ',diag_rdobs
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' LHN_RANGE = ',lhn_range
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' L_LHN = ',l_lhn
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,'(A,L1)')' L_LHN_1A = ',l_lhn_1a
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' L_LHN_SCALE = ',l_lhn_scale
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' L_LHN_SEARCH = ',l_lhn_search
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' LHN_DIAG = ',lhn_diag
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' EPSILON_LHN = ',epsilon_lhn
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' ALPHA_LHN = ',alpha_lhn
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' RELAX_CF_LHN = ',relax_cf_lhn
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' LHN_LIMIT = ',lhn_limit
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' L_VERIF_RANGE = ',l_verif_range
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' L_LHN_LIMIT = ',l_lhn_limit
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' L_LHN_FACT = ',l_lhn_fact
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' L_LHN_FILT = ',l_lhn_filt
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' FI_SCALE_LHN = ',fi_scale_lhn
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' NPASS_RF_LHN = ',npass_rf_lhn
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' REMOVE_NEG_LH = ',remove_neg_lh
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,*)' USE_CONV_IN_MOPS = ',use_conv_in_mops
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,'(A,I0)') ' a_assim_start_min = ',a_assim_start_min
CALL umPrint(lineBuffer,src='acp_namel')
WRITE(lineBuffer,'(A,I0)') ' a_assim_end_min = ',a_assim_end_min
CALL umPrint(lineBuffer,src='acp_namel')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='acp_namel')

END SUBROUTINE print_nlist_acp

! ------------------------------------------------------------------

SUBROUTINE read_nml_acp(nml_unit)

USE um_parcore, ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist, ONLY: setup_nml_type
USE parkind1, ONLY: jprb,jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Argument
INTEGER :: nml_unit

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_ACP'

INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 7 + modeacp + 6*nobtypmx
INTEGER, PARAMETER :: n_real = 5 + nobtypmx
INTEGER, PARAMETER :: n_log = 11

TYPE my_namelist
  SEQUENCE
  INTEGER :: ac_order     (nobtypmx)
  INTEGER :: ac_obs_types (nobtypmx)
  INTEGER :: no_obs_files
  INTEGER :: obs_format
  INTEGER :: lhn_range
  INTEGER :: npass_rf_lhn
  INTEGER :: diag_rdobs
  INTEGER :: macdiag(modeacp)
  INTEGER :: timea  (nobtypmx)
  INTEGER :: timeb  (nobtypmx)
  INTEGER :: tgetoba(nobtypmx)
  INTEGER :: tgetobb(nobtypmx)
  INTEGER :: a_assim_start_min
  INTEGER :: a_assim_end_min
  REAL :: epsilon_lhn
  REAL :: alpha_lhn
  REAL :: relax_cf_lhn
  REAL :: lhn_limit
  REAL :: fi_scale_lhn
  REAL :: nudge_lam (nobtypmx)
  LOGICAL :: lac_mes
  LOGICAL :: l_lhn
  LOGICAL :: l_lhn_1a
  LOGICAL :: l_lhn_scale
  LOGICAL :: l_lhn_search
  LOGICAL :: lhn_diag
  LOGICAL :: l_lhn_limit
  LOGICAL :: l_lhn_fact
  LOGICAL :: l_lhn_filt
  LOGICAL :: remove_neg_lh
  LOGICAL :: use_conv_in_mops
END TYPE my_namelist
        
TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,  &
                             n_real_in=n_real, n_log_in=n_log)

!     Namelist defaults for ACP

! Group dependent arrays in namelist.
ac_order(:)       = imdi
! Obs Type dependent arrays in namelist.
ac_obs_types(:)   = imdi
! Switch  for Mesoscale UM assimilation.
lac_mes = .FALSE.
!  Latent Heat Nudging
l_lhn = .FALSE.
! Use version 1A of  LHN
l_lhn_1a = .FALSE.
! Format of ac observation file
obs_format = imdi
! Number of observation files to be used
no_obs_files = imdi
! Group dependent arrays in namelist.
nudge_lam(:)      = rmdi
! Obs Type dependent arrays in namelist.
timeb(:)          = imdi
timea(:)          = imdi
tgetobb(:)        = imdi
tgetoba(:)        = imdi
!  Scale points by 1/EPSILON, in LHN, if no near rain found
l_lhn_scale = .FALSE.
!  Use search routine for nearby rain in LHN
l_lhn_search = .FALSE.
!  Display detailed diagnostics for LHN routine
lhn_diag = .FALSE.
!  Max range (in grid points) to search to in Latent Heat Nudging code
lhn_range = imdi
!  Limit size of increments
l_lhn_limit = .FALSE.
!  Maximum size of increment to theta within LHN (K/day)
lhn_limit = rmdi
!  Use limit in ALPHA_LHN
l_lhn_fact = .FALSE.
!  Filter LHN theta incrs
l_lhn_filt = .FALSE.
!  Minimum acceptable ratio of model ppn to observed, for LHN
epsilon_lhn = rmdi
!  Relaxation coefficient for theta increments from the LHN scheme
relax_cf_lhn = rmdi
!  Minimum factor by which model ppn can be scaled to fit obs ppn
alpha_lhn = rmdi
!  Recursive filter scale in metres
fi_scale_lhn = rmdi
!  Number of passes through recursive filter, for LHN increments
npass_rf_lhn = imdi
!  Remove negative LH
remove_neg_lh =  .FALSE.
!  Use convective precip
use_conv_in_mops = .FALSE.
! Control of diagnostic output from RDOBS/RDOBS2
! 0 - No listing ; 1 - Standard Listing ; 2 - More detailed listing.
diag_rdobs = imdi
! Set ac diagnostic mode (each timestep)
macdiag(:) = imdi

IF (mype == 0) THEN

  READ (UNIT=nml_unit, NML=acp, IOSTAT=errorstatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist ACP", iomessage)

  my_nml % ac_order(:)     = ac_order(:)
  my_nml % ac_obs_types(:) = ac_obs_types(:)
  my_nml % no_obs_files    = no_obs_files
  my_nml % obs_format      = obs_format
  my_nml % lhn_range       = lhn_range
  my_nml % npass_rf_lhn    = npass_rf_lhn
  my_nml % diag_rdobs      = diag_rdobs
  my_nml % macdiag(:)      = macdiag(:)
  my_nml % timea  (:)      = timea  (:)
  my_nml % timeb  (:)      = timeb  (:)
  my_nml % tgetoba(:)      = tgetoba(:)
  my_nml % tgetobb(:)      = tgetobb(:)
  my_nml % a_assim_start_min = a_assim_start_min
  my_nml % a_assim_end_min   = a_assim_end_min
  my_nml % epsilon_lhn     = epsilon_lhn
  my_nml % alpha_lhn       = alpha_lhn 
  my_nml % relax_cf_lhn    = relax_cf_lhn 
  my_nml % lhn_limit       = lhn_limit
  my_nml % fi_scale_lhn    = fi_scale_lhn
  my_nml % nudge_lam (:)   = nudge_lam (:)
  my_nml % lac_mes         = lac_mes
  my_nml % l_lhn           = l_lhn
  my_nml % l_lhn_1a        = l_lhn_1a
  my_nml % l_lhn_scale     = l_lhn_scale
  my_nml % l_lhn_search    = l_lhn_search
  my_nml % lhn_diag        = lhn_diag   
  my_nml % l_lhn_limit     = l_lhn_limit
  my_nml % l_lhn_fact      = l_lhn_fact
  my_nml % l_lhn_filt      = l_lhn_filt
  my_nml % remove_neg_lh   = remove_neg_lh  
  my_nml % use_conv_in_mops= use_conv_in_mops

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  ac_order(:)      = my_nml % ac_order(:)     
  ac_obs_types(:)  = my_nml % ac_obs_types(:)
  no_obs_files     = my_nml % no_obs_files
  obs_format       = my_nml % obs_format
  lhn_range        = my_nml % lhn_range
  npass_rf_lhn     = my_nml % npass_rf_lhn
  diag_rdobs       = my_nml % diag_rdobs
  macdiag(:)       = my_nml % macdiag(:)
  timea  (:)       = my_nml % timea  (:)
  timeb  (:)       = my_nml % timeb  (:)
  tgetoba(:)       = my_nml % tgetoba(:)
  tgetobb(:)       = my_nml % tgetobb(:)
  a_assim_start_min= my_nml % a_assim_start_min
  a_assim_end_min  = my_nml % a_assim_end_min
  epsilon_lhn      = my_nml % epsilon_lhn
  alpha_lhn        = my_nml % alpha_lhn
  relax_cf_lhn     = my_nml % relax_cf_lhn
  lhn_limit        = my_nml % lhn_limit
  fi_scale_lhn     = my_nml % fi_scale_lhn
  nudge_lam (:)    = my_nml % nudge_lam (:)
  lac_mes          = my_nml % lac_mes
  l_lhn            = my_nml % l_lhn
  l_lhn_1a         = my_nml % l_lhn_1a
  l_lhn_scale      = my_nml % l_lhn_scale
  l_lhn_search     = my_nml % l_lhn_search
  lhn_diag         = my_nml % lhn_diag
  l_lhn_limit      = my_nml % l_lhn_limit
  l_lhn_fact       = my_nml % l_lhn_fact
  l_lhn_filt       = my_nml % l_lhn_filt
  remove_neg_lh    = my_nml % remove_neg_lh
  use_conv_in_mops = my_nml % use_conv_in_mops

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_acp

SUBROUTINE check_nml_acp()

! Description:
!   Subroutine to check some control variables in AC scheme acp namelist.

USE chk_opts_mod, ONLY: chk_var, def_src
USE parkind1, ONLY: jprb,jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER (LEN=*),  PARAMETER :: RoutineName = 'CHECK_NML_ACP'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

IF (l_ac .AND. a_assim_start_min /= imdi) THEN
! Check sensible values chosen
  CALL chk_var( a_assim_start_min, 'a_assim_start_min', '[0:99999999]' )
  CALL chk_var( a_assim_end_min,   'a_assim_end_min',   '[0:99999999]' )
! These appear to be arrays but not clear whether
! more than 1 member actually used these days, nor how to check them
  CALL chk_var( ac_obs_types(1),   'ac_obs_types(1)',   '[1:9999]' )
  CALL chk_var( ac_order(1),       'ac_order(1)',       '[1:9999]' )
  CALL chk_var( alpha_lhn,         'alpha_lhn',         '[0.0:2.0]' )
  CALL chk_var( diag_rdobs,        'diag_rdobs',        '[0,1,2]' )
  CALL chk_var( epsilon_lhn,       'epsilon_lhn',       '[0.0:2.0]' )
  CALL chk_var( fi_scale_lhn,      'fi_scale_lhn',      '[0.0:100000.0]' )
  IF (l_lhn_limit) &  ! will be rmdi if l_lhn_limit == .false.
  CALL chk_var( lhn_limit,         'lhn_limit',         '[0.0:10.0]' )
  CALL chk_var( lhn_range,         'lhn_range',         '[0:1000]' )
  CALL chk_var( macdiag(1),        'macdiag(1)',        '[0:99]' )
  CALL chk_var( no_obs_files,      'no_obs_files',      '[0:10]' )
  CALL chk_var( npass_rf_lhn,      'npass_rf_lhn',      '[0:100]' )
  CALL chk_var( nudge_lam(1),      'nudge_lam(1)',      '[0.0:100000000.0]' )
  CALL chk_var( obs_format,        'obs_format',        '[2,3]' )
  CALL chk_var( relax_cf_lhn,      'relax_cf_lhn',      '[0.0:1.0]' )
  CALL chk_var( tgetoba(1),        'tgetoba(1)',        '[0:999999]' )
  CALL chk_var( tgetobb(1),        'tgetobb(1)',        '[0:999999]' )
  CALL chk_var( timea(1),          'timea(1)',          '[0:999999]' )
  CALL chk_var( timeb(1),          'timeb(1)',          '[0:999999]' )
END IF

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nml_acp

END MODULE acp_namel_mod
