PP field codes list
===================

The PP codes or field codes listed in this document are a code list for phenomenon identification which predates STASH.
In the STASHmaster file, the PP code is referred to as PPFC, in a PP or FieldsFile, the PP code is stored in LBFC.  In a PP or FieldsFile the valid values for the vertical coordinate type (LBVC) are also taken from this list.

.. Table start. keep this line just above the code table to help automatic parsing.

====  ==================================================================================
Code  Description
====  ==================================================================================
   0  Unspecified
   1  Height
   2  Depth
   3  Geopotential (= g*height)
   4  ICAO height
   5  Boundary layer height
   6  Non-dimensional soil model level
   7  Exner pressure
   8  Pressure
   9  Hybrid pressure
  10  Sigma (= p/surface p)
  11  T*u
  12  T*v
  13  T**2
  14  T*omega
  15  Height**2
  16  Temperature T
  17  Dew point temperature
  18  Dew point depression
  19  Potential temperature
  20  Maximum temperature
  21  Minimum temperature
  22  Wet bulb potential temperature
  23  Soil temperature (levels 2,3,4)
  24  d theta/dt
  25  Visibility
  26  Brunt-Vaisala Frequency 'N'
  27  (Atmospheric) density
  28  d(p*)/dt .  p* = surface pressure
  29  Cloud fraction below a level in feet
  30  Total cloud       See note a.
  31  High cloud        See note a.
  32  Medium cloud      See note a.
  33  Low cloud         See note a.
  34  Convective cloud  See note a.
  35  Contrail
  36  Fractional land amount
  37  Fractional sea-ice cover
  38  Atmospheric model land/sea mask   Land=1, Sea=0
  39  Coriolis parameter
  40  Omega  (= dp/dt)
  41  Wet bulb temperature
  42  Vertical velocity  (= dz/dt)
  43  Eta dot
  44  Time in seconds
  45  d(sigma)/dt
  46  u*q
  47  v*q
  48  'x' wind component (with respect to grid)
  49  'y' wind component (with respect to grid)
  50  Wind speed
  51  Thermal wind speed
  52  Vertical speed shear
  53  u*omega
  54  v*omega
  55  Wind direction
  56  Westerly component of wind  u
  57  Southerly component of wind v
  58  u**2
  59  v**2
  60  Product of wind components  u*v
  61  'x' component of wind stress
  62  'y' component of wind stress
  63  Kinetic energy
  64  Wind in 'DDDFFF' format
  65  Hybrid height
  66  Reserved for M.Mawson, RR will advise, 2000
  67  Reserved for M.Mawson, RR will advise, 2000
  68  u-acceleration from saturated stress
  69  v-acceleration from saturated stress
  70  u-acceleration from hydraulic jump
  71  v-acceleration from hydraulic jump
  72  Absolute vorticity
  73  Relative vorticity
  74  Divergence
  75  Height of base of lowest cloud in feet
  76  Height of top of lowest cloud in feet
  77  Total precipitation accumulation
  78  QCF - Frozen cloud water
  79  QCL - Liquid cloud water
  80  Stream function
  81  Velocity potential
  82  Ertel potential vorticity (Q)
  83  Thermal vorticity
  84  Quasi-geostrophic potential vorticity
  85  Montgomery stream function
  86  Geostrophic absolute vorticity
  87  dq/dt
  88  Relative humidity
  89  Precipitable water
  90  Total precipitation (= 94 + 102 + 103)
  91  Equivalent ice depth
  92  Actual ice depth
  93  Snow depth (water equivalent)
  94  Convective rain
  95  Specific humidity q
  96  Condensed water (per unit area)
  97  Total rainfall rate
  98  Convective rainfall rate
  99  Dynamic rainfall rate
 100  Local convective rainfall
 101  Mixing ratio
 102  Large scale rain
 103  Snow
 104  Total rain  (= 94 + 102)
 105  Evaporation
 106  Total soil moisture content (levels 1-4)
 107  Sublimation
 108  Snowfall rate mm/s
 109  Total runoff
 110  Snowmelt
 111  'Quick' runoff
 112  'Slow' runoff
 113  (Precipitation minus evaporation) rate (= 90 - 105)
 114  Evaporation rate
 115  Evaporation from soil surface
 116  Large scale snowfall
 117  Convective snowfall
 118  Large scale snowfall rate mm/s
 119  Convective snowfall rate mm/s
 120  Total large-scale precipitation (102+116)
 121  Total convective precipitation (94+117)
 122  Soil moisture in each (of 4) levels
 123  u-acceleration from trapped lee waves
 124  v-acceleration from trapped lee waves
 125  Vertical transmission coefficient
 126  Max C.A.T. level
 127  Sea bed level
 128  Mean sea level
 129  Surface
 130  Tropopause level
 131  Maximum wind level
 132  Freezing level
 133  Top of atmosphere
 134  -20 deg.C level 135  Upper level (height)
 136  Lower level (height)
 137  Upper level (pressure)
 138  Lower level (pressure)
 139  Wet bulb freezing level height (asl) m
 140  Salinity
 141  Snowmelt heat flux
 142  Upper level (hybrid)
 143  Lower level (hybrid)
 144  Unified model test diagnostic for output no. 1
 145  Unified model test diagnostic for output no. 2
 146  Unified model test diagnostic for output no. 3
 147  Unified model test diagnostic for output no. 4
 148  'X' component of geostrophic wind
 149  'Y' component of geostrophic wind
 150  Standard deviation of orography
 151  Distance to the centre of the earth
 152  Orography XX gradient component
 153  Orography XY gradient component
 154  Orography YY gradient component
 155  Orographic roughness
 156  Total ozone climatology
 157  Ultra violet index climatology
 158  Total ozone field, Dobson units
 159  UVB (Ultra violet index - dimensionless)
 160  Drag coefficient CD
 161  Theta_e
 162  Sunshine hours
 163  Convective cloud liquid re * convective cloud amount
 164  Layer cloud liquid re * layer cloud amount
 165  Convective cloud amount in SWRAD (microphysics)
 166  Layer cloud amount in SWRAD (microphysics)
 167  Layer cloud condensed water path * amount
 168  Cloud emissivity * cloud fraction
 169  Cloud albedo * cloud fraction
 170  Transmissivity
 171  RHOKH (RHO* * CH * SURF_LAYER_WIND_SHEAR)
 172  RHOKM (RHO* * CD * SURF_LAYER_WIND_SHEAR)
 173  Probability of visibility less than 5 km
 174  Silhouette area of unresolved orography per unit horizontal area
 175  Peak to trough height of unresolved orography divided by 2 (metres
 176  Latitude  (north positive)
 177  Longitude  (east positive)
 178  Sensible heat flux
 179  Soil heat flux
 180  Latent heat flux
 181  Bulk Richardson number
 182  Wind mixing energy
 183  CH - bulk transfer coefficient of heat
 184  Moisture flux
 185  Mass flux
 186  Net short wave radiation flux
 187  Net long wave radiation flux
 188  Total surface heat flux (inc. sensible+latent)
 189  Thermal advection
 190  C.A.T. probability
 191  Snow probability
 192  Boundary mixing coeffs.
 193  Convective heating rate
 194  Convective moistening rate
 195  Vertical momentum flux - U
 196  Vertical momentum flux - V
 197  Snow melt heating flux
 198  Evaporation duct height (previously called just 'Duct height')
 199  Evap. duct intensity (max wavelength trapped) (p.c.j. 'Duct intens
 200  Downward solar
 201  Upward solar
 202  Net surface radiation flux (No sensible or latent)
 203  Total downward surface solar flux over sea-ice.
 204  Reserved for M.Mawson, RR will advise, 2000
 205  IR down
 206  IR up
 207  Clear-sky flux (type II) solar up
 208  Clear-sky flux (type II) solar down
 209  Sea-ice temperature
 210  Clear-sky (type II) IR up
 211  Clear-sky (type II) IR down
 212  Convective cloud base times amount
 213  Dilute convectively available potential energy. J/kg (or equiv. m2
 214  Clear-sky (type II) net shortwave flux
 215  Clear-sky (type II) net longwave flux
 216  Cloud top height (asl)  kft
 217  Convective cloud top times amount
 218  Convective cloud water
 219  Convective cloud condensed water path kg/m**2
 220  Cloud amount in a layer / at a level  (reworded 4/3/93)
 221  Total water condensed by convection
 222  Convective cloud base level-number
 223  Convective cloud top level-number
 224  Advecting u * layer thickness (pascals) = zonal mass-flux
 225  Advecting v * layer thickness (pascals) = meridional mass-flux
 226  Blocking index (Tibaldi & Molteni definition)
 227  'EMAX' from convection scheme
 228  Total net radiative heating
 229  Zonal mass-flux * temperature
 230  Meridional mass-flux * temperature
 231  Zonal mass-flux * specific humidity q
 232  Meridional mass-flux * specific humidity q
 233  Zonal mass-flux * liquid water temperature
 234  Meridional mass-flux * liquid water temperature
 235  Zonal mass-flux * total water
 236  Meridional mass-flux * total water
 237  Clear-sky (type I) solar up
 238  Clear-sky (type I) solar down
 239  Vapour pressure
 240  Clear-sky (type I) IR up
 241  Clear-sky (type I) IR down
 242  Probability of precipitation
 243  Clear-sky (type I) net shortwave flux
 244  Clear-sky (type I) net longwave flux
 245  Clear-sky(type I) radiative heating
 246  Clear-sky(type I) short-wave radiative heating
 247  Clear-sky(type I) long wave radiative heating
 248  'X' component of ageostrophic wind
 249  'Y' component of ageostrophic wind
 250  Clear-sky(type II) radiative heating
 251  Short-wave radiative heating
 252  Clear-sky(type II) ditto
 253  Long-wave radiative heating
 254  Clear-sky(type II) ditto
 255  Convective cloud liquid re * convective cloud weighting (full leve
 256  Convective cloud weighting for microphysics (full levels)
 257  RHO_CD_MODV1 = rhostar*cD*modv1
 258  RHO_KM = rho*Km (Km = turbulent mixing coefficient for momentum)
 259  Atmospheric energy correction    W/m**2  (from UM 3.3)
 260  Sea-ice topmelt w/m2
 261  Sea-ice botmelt w/m2
 262  Fractional time of change of ice cover.
 263  Zonal mass-flux * u
 264  Meridional mass-flux * u
 265  Zonal mass-flux * v
 266  Meridional mass-flux * v
 267  Zonal mass-flux * geopotential
 268  Meridional mass-flux * geopotential
 269  Zonal mass-flux * moist static energy
 270  Meridional mass-flux * moist static energy
 271  Canopy water content
 272  Canopy condensation
 273  Canopy evaporation
 274  Canopy water throughfall
 275  Canopy height (all vegetation)
 276  Air concentration of radioactivity Becquerels(Bq) m**-3
 277  Dosage of radioactivity Bq.seconds m**-3
 278  Dry deposition of radioactivity Bq m**-2
 279  Wet deposition of radioactivity Bq m**-2
 280  Leads;net solar radiation wm-2
 281  Leads;net infra-red flux  wm-2
 282  Leads;sensible heat flux  wm-2
 283  Leads;latent heat flux    wm-2
 284  Total deposition of radioactivity Bq m**-2
 285  Evap from sea * leads fraction.
 286  Total aerosol (micro g/kg)
 287  Total aerosol emissions (micro g/sq m/s)
 288  Visibility assimilation weights
 289  Visibility assimilation increment
 290  P* (surface pressure) weights
 291  Theta weights
 292  Wind weights
 293  Surface wind weights.......
 294  RH weights
 295  Precipitation rate assimilation weights
 296  Spare for A/C
 297  Spare for A/C
 298  Spare for A/C
 299  Specific layer cloud water content
 300  Spare for A/C
 301  Spare for A/C
 302  Spare for A/C
 303  P* (surface pressure) assimilation increments
 304  Theta assimilation increments
 305  u assimilation increments
 306  v assimilation increments
 307  Hydrostatic increments
 308  Probability of ground frost
 309  RH assimilation increments
 310  Precipitation rate assimilation increments
 311  Reserved for M.Mawson, RR will advise, 2000
 312  Reserved for M.Mawson, RR will advise, 2000
 313  Reserved for M.Mawson, RR will advise, 2000
 314  Reserved for M.Mawson, RR will advise, 2000
 315  Increment in q from a routine
 316  Increment in T from a routine
 317  Increment in theta from a routine
 318  Increment in u from a routine
 319  Increment in v from a routine
 320  Surface roughness (heat)
 321  Root depth
 322  Snow free albedo (all types of land)
 323  Surface resistance to evaporation
 324  Surface roughness (momentum)
 325  Surface capacity
 326  Vegetation fraction
 327  Veg. infilt. enhancement factor
 328  Deep snow albedo
 329  Wilting point
 330  Critical point
 331  Field capacity
 332  Saturation
 333  Saturated conductivity
 334  Eagleson's exponent
 335  Heat capacity
 336  Heat conductivity
 337  SOILB (Soil hydrology parameter 'BS')
 338  Field capacity-root depth
 339  Surface infilt.-infilt. factor
 340  Thermal inertia-(lambda C ** 1/2)
 341  Veg. infilt. enhancement factor * Saturated conductivity (327*333)
 342  Saturated soil water suction (m)
 343  Lowest convective cloud amount
 344  Lowest convective cloud base (pa)
 345  Lowest convective cloud top (pa)
 346  Lowest convective cloud base (Kft)
 347  Lowest convective cloud top (Kft)
 348  'X' component of 'Q'
 349  'Y' component of 'Q'
 350  DIV(Q)
 351  Wave/ Spectral wave energy component
 352  Wave/
 353  Wave/ Energy increments after propagation
 354  Wave/ Energy increments after input source term
 355  Wave/ Energy increments after nonlinear transfer
 356  Wave/ Energy increments after dissipation
 357  Wave/ Energy increments after bottom friction
 358  Wave/ Energy increments after assimilation
 359  Wave/ Assimilation increments for wave height
 360  Wave/ Assimilation increments for windspeed
 364  Wave/ Wave induced surface stress x component
 365  Wave/ Wave induced surface stress y component
 366  Wave/ Induced stress
 367  Wave/ Dependent roughness length
 368  Wave/ Total induced current speed
 369  Wave/ Total induced current direction
 370  Wave/ Surge induced current speed
 371  Wave/ Surge induced current direction
 372  Wave/ Tide induced water level
 373  Wave/ Tide induced current speed
 374  Wave/ Tide induced current direction
 375  Wave/ Mean wave height error
 376  Wave/ Mean wave period error
 377  Wave/ Mean wave speed error
 378  Wave/ RMS wave height error
 379  Wave/ RMS wave period error
 380  Wave/ RMS wave speed error
 384  Wave/ Water temperature
 385  Wave/ Height of wind-driven waves
 386  Wave/ Swell height
 387  Wave/ Combined wave and swell height
 388  Wave/ Windsea upcrossing period
 389  Wave/ Windsea principal direction
 390  Wave/ Swell upcrossing period
 391  Wave/ Swell principal direction
 392  Wave/ Total depth (wave+surge)
 393  Wave/ Total wave upcrossing period
 394  Wave/ Total wave principal direction
 395  Wave/ Grid point type (sea/land/coast)
 396  Wave/ Surge induced water level
 397  Wave/ Wave model land/sea indicator
 401  MA/ d(ln p)/dt
 402  MA/ d(theta)/dp
 403  MA/ Log(Absolute ertel pot. vorticity)
 404  MA/ Sign of dq/dt
 405  MA/ d(theta)/dt
 406  MA/ Shortwave heating rates
 407  MA/ Longwave cooling rates
 408  MA/ Total radiative heating
 409  MA/ Friction term in circulation budget
 410  MA/ dp/d(theta)
 411  MA/ Geostrophic Ertel potential vorticity
 412  MA/ dq/dt advection term
 413  MA/ dq/dt source term
 414  MA/ Total dq/dt ie source + advection
 415  MA/ d(theta)/dt
 416  MA/ Absolute vorticity (interpolated)
 417  MA/ Absolute vorticity*total dq/dt
 418  MA/ dq/dt from Raleigh friction
 420  EASTWARD FLUX - SPECTRAL PSEUDOMOM
 421  SOUTHWARD FLUX - SPECTRAL PSEUDOMOM
 422  WESTWARD FLUX - SPECTRAL PSEUDOMOM
 423  NORTHWARD FLUX - SPECTRAL PSEUDOMOM
 424  EASTWARD FORCE FROM SPECTRAL GW
 425  NORTHWARD FORCE FROM SPECTRAL GW
 451  MAtr/ Passive tracer
 452  MAtr/ Nitrous oxide tracer
 453  MAtr/ Ozone tracer
 501  MAUA/ Atmosphere tracer 1 (conventionally O3)
 502  MAUA/ Atmosphere tracer 2 (conventionally H2O)
 503  MAUA/ Atmosphere tracer 3 (conventionally CO)
 504  MAUA/ Atmosphere tracer 4 (conventionally CH4)
 505  MAUA/ Atmosphere tracer 5 (conventionally N2O)
 506  MAUA/ Atmosphere tracer 6 (conventionally NO)
 507  MAUA/ Atmosphere tracer 7 (conventionally NO2)
 508  MAUA/ Atmosphere tracer 8 (conventionally HNO3)
 509  MAUA/ Atmosphere tracer 9 (conventionally N2O5)
 510  MAUA/ Atmosphere tracer 10 (conventionally ClONO2)
 511  MAUA/ Atmosphere tracer 11 (conventionally ClO)
 512  MAUA/ Atmosphere tracer 12 (conventionally HCl)
 513  MAUA/ Atmosphere tracer 13 (conventionally CF2Cl2)
 514  MAUA/ Atmosphere tracer 14 (conventionally CFCl3)
 515  MAUA/ Atmosphere tracer 15 (conventionally HF)
 516  MAUA/ Atmosphere tracer 16
 517  MAUA/ Atmosphere tracer 17
 518  MAUA/ Atmosphere tracer 18
 519  MAUA/ Atmosphere tracer 19 (conventionally SO2)
 520  MAUA/ Atmosphere tracer 20 (conventionally DMS)
 521  MAUA/ Atmosphere tracer 21 (conventionally H2S)
 522  MAUA/ Atmosphere tracer 22 (conventionally WATER-SOLUBLE)
 523  MAUA/ Atmosphere tracer 23 (conventionally DUST-LIKE)
 524  MAUA/ Atmosphere tracer 24 (conventionally OCEANIC)
 525  MAUA/ Atmosphere tracer 25 (conventionally SOOT)
 526  MAUA/ Atmosphere tracer 26 (conventionally VOLCANIC ASH)
 527  MAUA/ Atmosphere tracer 27 (conventionally SULPHURIC ACID)
 528  MAUA/ Atmosphere tracer 28 (conventionally AMMONIUM SULPHATE)
 529  MAUA/ Atmosphere tracer 29 (conventionally MINERAL)
 530  MAUA/ Atmosphere tracer 30 (conventionally OCFF)
 568  MAUA/
 569  MAUA/ Sulphur Dioxide emissions.
 570  MAUA/ Dimethyl Sulphide emissions.
 571  MAUA/ Hydrogen Sulphide emissions.
 572  MAUA/ NH3 gas emissions kg/m2/s
 573  MAUA/ Soot emissions   kg/m2/s
 574  MAUA/ Surface emissions of biomass smoke
 575  MAUA/ Oraganic carbon from fossil fuels (OCFF)
 576  DIMETHYL SULPHIDE EMISSIONS kg(S) m-2 s-1
 580  MAUA/ Hydroxyl radical concentration molecules cm**-3
 581  MAUA/ Hydrogen Peroxide concentration
 600  MAUA/ HO2 concentration
 601  OCtr/ Temperature.
 602  OCtr/ Salinity.
 603  OCtr/ T*
 604  OCtr/ Temperature at t-1.
 605  OCtr/ Salinity at t-1.
 606  OCtr/ No. deep convection points
 607  OCtr/ Fractional expansion of sea water (dimensionless)
 608  OCtr/ Sea surface elevation
 609  OCtr/ Density of sea water/kg m**-3
 610  OCtr/ Sea level rise/m
 611  OCtr/ Stream function.
 612  OCtr/ Stream function at t-1.
 613  OCtr/ Change of stream function across a time step.
 614  OCtr/ Change of stream function across previous timestep.
 615  OCtr/ (Formerly reciprocal of total depth at U/V points, use 715 1
 616  OCtr/ Number of vertical levels of ocean at T points.
 617  Octr/ rigid-lid pressure/Pa
 618  OCtr/ Change of vorticity across one timestep. Units s*-2 from Jun
 619  OCtr/ Vertical expansion/m of sea water.
 620  OCtr/ Convergence heat WM-2
 621  OCtr/ Equivalent ice depth.
 622  OCtr/ Snow depth.
 623  OCtr/ Snowfall.
 624  OCtr/ Sublimation.
 625  OCtr/ Net solar heat flux.
 626  OCtr/ Net surface heat flux ('HTN').
 627  OCtr/ Wind mixing energy.
 628  OCtr/ Diffusive heat flux.
 629  OCtr/ Precipitation minus evaporation ('PLE')
 630  OCtr/ Sum of net solar and net surface heating.
 631  OCtr/ River outflow
 632  OCtr/ Surface water flux * salinity / density  m s**-1
 633  OCtr/ 'x' component of sea-ice/ocean surface stress ('ISX')
 634  OCtr/ 'y' component of sea-ice/ocean surface stress ('ISY')
 635  OCtr/ Buoyancy/kg m**-2 s**-2
 636  OCtr/ Buoyancy flux/kg m**-1 s**-3
 637  OCtr/ In-situ temperature
 638  OCtr/ sigma-T/kg m**-3 ( (rho(insitu_T,S,0)-1000 )
 639  OCtr/ sigma-theta/kg m**-3 ( rho(theta,S,0)-1000 )
 640  OCtr/ Zonal heat advection.
 641  OCtr/ Meridional heat advection.
 642  OCtr/ Heating rate due to advection  K s**-1
 643  OCtr/ Heating rate due to diffusion  K s**-1
 644  OCtr/ Heating rate due to surface fluxes  K s**-1
 645  OCtr/ Heating rate due to mixing  K s**-1
 646  OCtr/ Heating rate due to filtering  K s**-1
 647  OCtr/ Heating rate  K s**-1
 648  OCtr/ Rate of change of salinity/s**-1 (UM can't control LBPROC)
 649  OCtr/ Climatological reference surface salinity
 650  OCtr/ Climatological reference SST (sea surface temperature).
 651  OCtr/ Anomalous heat flux from Haney term.
 652  OCtr/ Anomalous latent heat flux from Haney term
 653  OCtr/ Mixed layer depth,type 1.
 654  OCtr/ Mixed layer depth,type 2.
 655  OCtr/ Heat content /J
 656  OCtr/ Latent heat content /J
 657  OCtr/ Turbulent kinetic energy m**2 s**-2
 658  OCtr/ Thickness diff coeff (ocean) cm2/s
 659  OCtr/
 660  OCtr/ Vertical mean vorticity forcing: advection s*-2
 661  OCtr/ Vertical mean vorticity forcing: hor diffn s*-2
 662  OCtr/ Vertical mean vorticity forcing: vrt diffn s*-2
 663  OCtr/ Vertical mean vorticity forcing: coriolis s*-2
 664  OCtr/ Vertical mean vorticity forcing: pressure s*-2
 665  OCtr/ Vertical integral vorticity forcing: advection cm s*-2
 666  OCtr/ Vertical integral vorticity forcing: hor diffn cm s*-2
 667  OCtr/ Vertical integral vorticity forcing: vrt diffn cm s*-2
 668  OCtr/ Vertical integral vorticity forcing: coriolis cm s*-2
 669  OCtr/ Vertical integral vorticity forcing: bottom p cm s*-2
 670  OCtr/ Boundary profiles.
 671  OCtr/ Anomalous heat flux.
 672  OCtr/ Anomalous salt flux.
 673  OCtr/ Virtual air-sea flux of co2
 674  OCtr/ Virtual air-sea flux of alkalinity
 675  OCtr/ Climatological Ice depth
 676  OCtr/ Mass of ocean below rigid lid/kg m**-2
 677  OCtr/ Ocean bottom pressure/Pa
 678  OCtr/ Anomolous sea-ice heat flux W,M2
 679  OCtr/ Anomalous sea-ice P-E flux kg, M2, S
 680  OCtr/ W (vertical velocity).
 681  OCtr/ Topmelt.
 682  OCtr/ Botmelt.
 683  OCtr/ Fractional ice cover.
 684  OCtr/ Ocean-ice heat flux.
 685  OCtr/ Carry heat.
 686  OCtr/ Carry salt.
 687  OCtr/ Actual ice depth.
 688  OCtr/ Snow depth over ice (m).
 689  OCtr/ Sea-ice internal pressure (usually in N.m-2)
 690  OCtr/ Sea-ice strength     (usually in N.m-2)
 692  OCtr/ Average snow depth.
 693  OCtr/ Average ice cover.
 694  OCtr/ Average O/I heat flux.
 695  OCtr/ KT Vertical tracer diffusion coeff
 696  OCtr/ RiT Richardson no tracer
 697  OCtr/ HT Max depth Large scheme tracer
 698  OCtr/ Downwards solar radiation over ice
 699  OCtr/ Depth/m of ocean at TS points
 700  OCuv/ River run-off ocean entry point co-ordinates.
 701  OCuv/ Baroclinic component of zonal velocity
 702  OCuv/ Baroclinic component of meridional velocity
 703  OCuv/ Zonal velocity (total)
 704  OCuv/ Meridional velocity (total)
 711  OCuv/ Barotropic component of zonal velocity
 712  OCuv/ Barotropic component of meridional velocity
 713  OCuv/ Baroclinic x-acceleration cm s**-2
 714  OCuv/ Baroclinic y-acceleration cm s**-2
 715  OCuv/ Reciprocal of total depth
 717  OCuv/ Number of vertical levels of ocean at U/V points.
 718  OCuv/ Vertical integral of zonal velocity/m2 s-1
 719  OCuv/ Vertical integral of meridional velocity/m2 s-1
 721  OCuv/ Zonal component of windstress ('TAUX')
 722  OCuv/ Meridional component of windstress ('TAUY')
 728  OCuv/ 'x' component of sea-ice velocity
 729  OCuv/ 'y' component of sea-ice velocity
 730  OCuv/ Magnitude of sea-ice velocity
 731  OCuv/ Zonal component of sea-ice velocity
 732  OCuv/ Meridional component of sea-ice velocity
 733  OCuv/ Zonal component of ice-ocean stress/Pa
 734  OCuv/ Meridional component of ice-ocean stress/Pa
 735  OCuv/ Zonal component of Coriolis stress/Pa
 736  OCuv/ Meridional component of Coriolis stress/Pa
 737  OCuv/ Zonal component of internal ice stress /Pa
 738  OCuv/ Meridional component of internal ice stress /Pa
 740  OCuv/ Zonal mean tracer transport diagnostics
 741  OCuv/ Baroclinic acceleration diag, meridional: Advection
 742  OCuv/ Baroclinic acceleration diag, meridional: Hor Diffusion
 743  OCuv/ Baroclinic acceleration diag, meridional: Vert Diffusion
 744  OCuv/ Baroclinic acceleration diag, meridional: Coriolis
 745  OCuv/ Baroclinic acceleration diag, meridional: Pressure
 746  OCuv/ Baroclinic acceleration diag, zonal: Advection
 747  OCuv/ Baroclinic acceleration diag, zonal: Hor Diffusion
 748  OCuv/ Baroclinic acceleration diag, zonal: Vert Diffusion
 749  OCuv/ Baroclinic acceleration diag, zonal: Coriolis
 750  OCuv/ Baroclinic acceleration diag, zonal: Pressure
 751  OCuv/ Meridional net pressure gradient
 752  OCuv/ Zonal net pressure gradient
 753  OCuv/ Meridional acceleration
 754  OCuv/ Zonal acceleration
 795  OCuv/ KM Vertical momentum diffusion coeff
 796  OCuv/ RiM Richardson no momentum
 797  OCuv/ HM Max depth Large scheme momentum
 801  OCE/ ocean extra tracer 1 (conventionally TCO2 )
 802  OCE/ ocean extra tracer 2 (conventionally alkalinity)
 803  OCE/ ocean extra tracer 3 (conventionally nutrient)
 804  OCE/ ocean extra tracer 4 (conventionally phytoplankton conc.)
 805  OCE/ ocean extra tracer 5 (conventionally zooplankton)
 806  OCE/ ocean extra tracer 6 (conventionally detritus)
 807  OCE/ ocean extra tracer 7 (conventionally tritium)
 808  OCE/ ocean extra tracer 8 (conventionally 3H+3He total mass)
 809  OCE/ ocean extra tracer 9 (conventionally CFC11)
 810  OCE/ ocean extra tracer 10 (conventionally CFC12)
 811  OCE/ ocean extra tracer 11 (conventionally CFC13)
 812  OCE/ ocean extra tracer 12 (conventionally carbon 14)
 813  OCE/ ocean extra tracer 13
 814  OCE/ ocean extra tracer 14
 815  OCE/ ocean extra tracer 15
 816  OCE/ ocean extra tracer 16
 817  OCE/ ocean extra tracer 17
 818  OCE/ ocean extra tracer 18
 819  OCE/
 820  OCE/ Sonic layer depth (m)
 821  OCE/ Sonic layer trap at 10Hz (0/1)
 822  OCE/ Sonic layer trap at 50Hz (0/1)
 823  OCE/ Sonic layer trap at 1kHz (0/1)
 824  OCE/ Sonic layer trap at 10kHz (0/1)
 825  OCE/ Sound channel strength at 10Hz (1-3)
 826  OCE/ Sound channel strength at 50Hz (1-3)
 827  OCE/ Sound channel strength at 1kHz (1-3)
 828  OCE/ Sound channel strength at 10kHz (1-3)
 829  OCE/ Range to the 1st CZ (km)
 830  OCE/ Strength of 1st CZ (0-100)
 831  OCE/ Gent and McWilliams scheme eddy u velocity (ocean)
 832  OCE/ Gent and McWilliams scheme eddy v velocity (N face) (ocean)
 833  OCE/ Gent and McWilliams scheme eddy w velocity (top face) (ocean
 834  OCE/ d theta/dt from Gent and McWilliams scheme K/S
 835  OCE/ Speed of sound in water m/s
 836  OCE/ Depth of sound speed minimum (sound channel) m
 837  OCE/ Depth of max neg ssp grad (m)
 838  OCE/ CO2 atmosphere level 1 conc ppmv
 840  OCE/ Total temperature advection, zonal.
 841  OCE/ Total temperature diffusion, zonal.
 842  OCE/ Small sediment concentration (g m**-3)
 843  OCE/ Large sediment concentration (g m**-3)
 844  OCE/ Total sediment concentration (g m**-3)
 850  OA/ Surface height analysis weights
 851  OA/ Mixed layer depth analysis weights
 852  OA/ Surface temperature analysis weights
 853  OA/ Potential temperature analysis weights
 854  OA/ Salinity analysis weights
 855  OA/ Velocity components analysis weights
 860  OA/ Surface height analysis increments
 861  OA/ Mixed layer depth analysis increments
 862  OA/ Surface temperature analysis increments
 863  OA/ Potential temperature analysis increments
 864  OA/ Salinity analysis increments
 865  OA/ Zonal velocity analysis increments
 870  OA/ Meridional velocity analysis increments
 871  OA/ Meridional velocity increments after surface height analysis
 872  OA/ Meridional velocity increments after thermal analysis
 873  OA/ Meridional velocity increments after saline analysis
 876  OA/ Zonal velocity increments after surface height analysis
 877  OA/ Zonal velocity increments after thermal analysis
 878  OA/ Zonal velocity increments after saline analysis
 880  OA/ Pressure increments after surface height analysis
 881  OA/ Pressure increments after thermal analysis
 882  OA/ Pressure increments after saline analysis
 885  OA/ Potential temperature increments after surface height analys
 888  OA/ Saline increments after surface height analysis
 891  OCE/ PRIMARY PRODUCTION (GC/M2/DAY)
 892  OCE/ ZOOPLTN PRODUCTION (GC/M2/DAY)
 893  OCE/ PHYTO SPECIFIC GROWTH RATE (1/DAY)
 894  OCE/ PHYTO SPECIFIC GRAZING RATE (1/DAY)
 895  OCE/ PHYTO SPECIFIC MORTALITY (1/DAY)
 896  OCE/ NITRATE GAIN-EXCRETION (MMOL-N/M2/D)
 897  OCE/ NITRATE LOSS - GROWTH (MMOL-N/M2/D)
 898  OCE/ NITRATE GAIN-PHY MORT (MMOL-N/M2/D)
 899  OCE/ NITRATE GAIN-ZOO MORT (MMOL-N/M2/D)
 900  OCE/ NITRATE GAIN-PHY RESP (MMOL-N/M2/D)
 901  OCE/ NITRATE GAIN-REMIN    (MMOL-N/M2/D)
 902  OCE/ NUTRIENT LIMITATION
 903  OCE/ LIGHT LIMITATION
 904  OCE/ TEMPERATURE LIMITATION
 905  OCE/ DETRITUS FLUX  (MMOL-N/M2/D)
 906  OCE/ VERTICAL NITRATE FLUX  (MMOL-N/M2/D)
 907  OCE/ HORIZ NITRATE ADVECT RATE(MMOL/M3/D)
 908  OCE/ VERT NITRATE ADVECTN RATE(MMOL/M3/D)
 909  OCE/ HORIZ NITRATE DIFFUSION  (MMOL/M3/D)
 910  OCE/ VERT NITRATE DIFFUSION   (MMOL/M3/D)
 911  OCE/ NITRATE MIXING DUE TO MLM(MMOL/M3/D)
 912  OCE/ NITRATE CONVECTION       (MMOL/M3/D)
 913  OCE/ NITRATE CHANGE - BIOLOGY (MMOL/M3/D)
 914  OCE/ NITRATE CHANGE-RESETTING (MMOL/M3/D)
 915  OCE/ HORIZ PHYTO ADVECT RATE(MMOL-N/M3/D)
 916  OCE/ HORIZ ZOO   ADVECT RATE(MMOL-N/M3/D)
 917  OCE/ HORIZ DETRI ADVECT RATE(MMOL-N/M3/D)
 918  OCE/ Rate of change of sea-ice concentration /s-1
 919  OCE/ Rate of change of sea-ice depth /m s-1
 920  OCE/ Rate of change of snowdepth on sea-ice /m s-1
 921  OCE/Ocean Nr Surface Chlorophyll content (kg m-3)
 922  OCE/Freshwater flux from Iceberg Calving
 940  OCE/ Total temperature advection, meridional.
 941  OCE/ Total temperature diffusion, meridional.
 942  OCE/ Mmeridional overturning streamfunction/Sv
 943  Total freshwater flux into ocean at t-1 (kg/m**3)(cm/s)
 944  Unsmoothed surface pressure tendency (g/cm/s**2)
 945  Maximum current in water column (zonal component)
 946  Maximum current in water column (meridional component)
 947  Depth of maximum current in water column
1001  wqt flux
1002  wql flux
1003  wthetal flux
1004  wthetav flux
1005  sub cloud layer convective velocity scale
1006  cumulus layer convective velocity scale
1007  cloud base mass flux
1008  congestus indicator 1
1009  congestus indicator 2
1010  termination model level for congestus
1011  height of top of shallow convection
1012  height of base of shallow convection
1013  height of top of congestus convection
1014  height of base of congestus convection
1015  CONVECTIVE BOUNDARY LAYER VELOCITY  SCALE
1016  SURFACE BUOYANCY FLUX
1017  GRADIENT RICHARDSON NUMBER
1018  VERTICAL BUOYANCY GRADIENT
1019  MODULUS OF WIND SHEAR
1020  BL MOMENTUM DIFFUSION
1021  BL HEAT DIFFUSION
1022  TURBULENT KINETIC ENERGY
1023  X-COMP OF DIST OROGRAPHIC STRESS
1024  Y-COMP OF DIST OROGRAPHIC STRESS
1025  X-COMP SURFACE BL STRESS
1026  Y-COMP SURFACE BL STRESS
1027  Direct UV Flux (on 38 levels)
1028  Upward UV Flux (on 38 levels)
1029  Net Downward UV FLUX (on 38 levels)
1030  MEAN SOLAR BEARING OVER RAD TS
1031  SLOPE ASPECT
1032  SLOPE ANGLE
1033  OROG CORR FACTOR TO DIRECT SURF SW
1034  EXTRA SW SURF FLUX AFTER OROG CORR
1035  Mass flux on half levels (i.e. rho levels).
1036  COMBINED BOUNDARY LAYER TYPE DIAGNOSTIC
1037  nitrate (mmol m-3)
1038  phosphate (mmol m-3)
1039  diatom biomass (mg-C m-3)
1040  flagellate biomass (mg-C m-3)
1041  picoplankton biomass (mg-C m-3)
1042  dinoflagellate biomass (mg-C m-3)
1043  dissolved oxygen concentration (mmol m-3)
1044  zooplankton biomass (mg-C m-3)
1045  underwater visibility (m)
1046  daily averaged visibility (m)
1047  0/351 Clim Biogenic Aerosol
1048  0/352 Clim Biomass-burning (fresh)
1049  0/353 Clim Biomass-burning (aged)
1050  0/354 Clim Biomass-burning (in-cloud)
1051  0/355 Clim Black Carbon (fresh)
1052  0/356 Clim Black Carbon (aged)
1053  0/357 Clim Sea-salt (film mode)
1054  0/358 Clim Sea-salt (jet mode)
1055  0/359 Clim Sulphate (accumulation mode)
1056  0/360 Clim Sulphate (Aitken mode)
1057  0/361 Clim Sulphate (dissolved)
1058  0/362 Clim Dust size division 1
1059  0/363 Clim Dust size division 2
1060  0/364 Clim Dust size division 3
1061  0/365 Clim Dust size division 4
1062  0/366 Clim Dust size division 5
1063  0/367 Clim Dust size division 6
1064  0/368 Reserved for other aerosol clims
1065  0/369 Reserved for other aerosol clims
1066  0/370 Reserved for other aerosol clims
1067  0/371 Reserved for other aerosol clims
1068  0/372 Reserved for other aerosol clims
1069  0/373 Reserved for other aerosol clims
1070  0/374 Reserved for other aerosol clims
1071  2/289 Biogenic aerosol optical depth
1072  Reserved just in case
1073  Reserved just in case
1074  Grid-box mean roughness length for momentum without orographic enhancement (m).
1075  RESIDUAL MN MERID. CIRC. VSTARBAR
1076  RESIDUAL MN MERID. CIRC. WSTARBAR
1077  ELIASSEN-PALM FLUX (MERID. COMPNT)
1078  ELIASSEN-PALM FLUX (VERT. COMPNT)
1079  DIVERGENCE OF ELIASSEN-PALM FLUX
1080  Main ozone tracer:  prognostic ozone tracer
1081  Coeff 1:            car o3 prod-loss (P-L)
1082  Coeff 2:            car o3 P-L wrt o3 mix ratio
1083  Coeff 3:            car o3 vol mixing ratio
1084  Coeff 4:            car o3 P-L wrt temp
1085  Coeff 5:            car o3 clim temp
1086  Coeff 6:            car o3 P-L wrt o3 above point
1087  Coeff 7:            car o3 column above point
1088  reserved for C Mathison
1089  reserved for C Mathison
1090  Very Low Cloud amount
1091  Height of decoupled layer base (m)
1092  Height of stratocumulus cloud base (m)
1093  Parametrized entrainment rate for surface-based mixed layer (m/s)
1094  Parametrized entrainment rate at the boundary layer top (m/s)
1095  East. FLUX SPECTRAL PSEUDOMOM. P LEVS
1096  West. FLUX SPECTRAL PSEUDOMOM. P LEVS
1097  EAST. FORCE FROM SPECTRAL GW P LEVS
1098  X COMPT OF GRAV. WAVE STRESS P LEVS
1099  U-ACCEL FROM SATURATED STRESS P LEVS
1100  CCRad : Lowest conv. cloud base layer
1101  CCRad : CCW passed to radiation
1102  Surface emissivity at 23.8 GHz (F. Karbou)
1103  Standard deviation of surface emissivity at 23.8 GHz (F.Karbou)
1104  Surface emissivity at 31.4 GHz (F. Karbou)
1105  Standard deviation of surface emissivity at 31.4 GHz (F. Karbou)
1106  Surface emissivity at 50.0 GHz (F. Karbou)
1107  Standard deviation of surface emissivity at 50.0 GHz (F. Karbou)
1108  Surface emissivity at 89.0 GHz (F. Karbou)
1109  Standard deviation of surface emissivity at 89.0 GHz (F. Karbou)
1110  Surface emissivity at 150.0 GHz (F. Karbou)
1111  Standard deviation of surface emissivity at 150.0 GHz (F. Karbou)
1112  Shear vorticity
1113  Curvature vorticity
1114  Zenith Total delay (ZTD)
1115  Total wave directional spread
1116  Secondary swell field sig. height
1117  Tertiary swell field sig. height
1118  Secondary swell field period
1119  Tertiary swell field period
1120  Wind sea wavelength
1121  Primary swell field wavelength
1122  Secondary swell field wavelength
1123  Tertiary swell field wavelength
1124  Secondary swell field principle dir
1125  Tertiary swell field principle dir
1126  Wind sea directional spread
1127  Primary swell field dir. spread
1128  Secondary swell field dir. spread
1129  Tertiary swell field dir. spread
1130  Primary swell wind sea fraction
1131  Secondary swell wind sea fraction
1132  Teriary swell wind sea fraction
1133  Total wind sea fraction
1134  CAPE timescale,deep (secs).
1135  Indicator of reduced CAPE timescale.
1136  Distance from Coast (m)
1137  Radiative Screen Temp on Tiles (K)
1138  Radiative Screen Temp on Sea Ice (K)
1139  Time since Transition (s)
1140  ......  Free for General use  ......
1141  ......  Free for General use  ......
1142  ......  Free for General use  ......
1143  ......  Free for General use  ......
1144  ......  Free for General use  ......
1145  ......  Free for General use  ......
1146  ......  Free for General use  ......
1147  ......  Free for General use  ......
1148  ......  Free for General use  ......
1149  BL flux of atmospheric tracer 30 (conventionally OCFF)
1150  Mass mixing ratio fresh OCFF
1151  Mass mixing ratio aged OCFF
1152  Mass mixing ratio cloud OCFF
1153  Fresh OCFF surface dry deposition flux kg/m2/s
1154  Aged OCFF surface dry deposition flux kg/m2/s
1155  OCFF in cloud occult deposition flux kg/m2/s
1156  Large scale rainout of OCFF
1157  Large scale washout of OCFF
1158  Convective scavenging of OCFF
1160  RHOKH_MIX
1161  RHO_ARESIST (RHOSTAR*CD_STD*VSHR)
1162  ARESIST [ 1/(CD_STD*VSHR) ]
1163  RESIST_B (1/CH-1/CD_STD)/VSHR
1164  DTRDZ_CHARNEY_GRID
1165  GRID-LEVEL OF SML INVERSION
1166  Rho * entrainment rate
1167  Fraction of the timestep
1168  zrzi
1169  GRID-LEVEL OF DSC INVERSION
1170  Rho * entrainment rate dsc
1171  Fraction of the timestep dsc
1172  zrzi dsc
1173  ZHSC  Top of decoupled layer
1174  Surface layer resist for dust div1
1175  Surface layer resist for dust div2
1176  Surface layer resist for dust div3
1177  Surface layer resist for dust div4
1178  Surface layer resist for dust div5
1179  Surface layer resist for dust div6
1301  BL flux of atmospheric tracer  1 (conventionally O3)
1302  BL flux of atmospheric tracer  2 (conventionally H2O)
1303  BL flux of atmospheric tracer  3 (conventionally CO)
1304  BL flux of atmospheric tracer  4 (conventionally CH4)
1305  BL flux of atmospheric tracer  5 (conventionally N2O)
1306  BL flux of atmospheric tracer  6 (conventionally NO)
1307  BL flux of atmospheric tracer  7 (conventionally NO2)
1308  BL flux of atmospheric tracer  8 (conventionally HNO3)
1309  BL flux of atmospheric tracer  9 (conventionally N2O5)
1310  BL flux of atmospheric tracer 10 (conventionally ClONO2)
1311  BL flux of atmospheric tracer 11 (conventionally ClO)
1312  BL flux of atmospheric tracer 12 (conventionally HCl)
1313  BL flux of atmospheric tracer 13 (conventionally CF2Cl2)
1314  BL flux of atmospheric tracer 14 (conventionally CFCl3)
1315  BL flux of atmospheric tracer 15 (conventionally HF)
1316  BL flux of atmospheric tracer 16
1317  BL flux of atmospheric tracer 17
1318  BL flux of atmospheric tracer 18
1319  BL flux of atmospheric tracer 19 (conventionally SO2)
1320  BL flux of atmospheric tracer 20 (conventionally DMS)
1321  BL flux of atmospheric tracer 20 (conventionally DMS)
1322  BL flux of atmospheric tracer 21 (conventionally H2S)
1323  BL flux of atmospheric tracer 22 (conventionally WATER-SOLUBLE)
1324  BL flux of atmospheric tracer 23 (conventionally DUST-LIKE)
1325  BL flux of atmospheric tracer 24 (conventionally OCEANIC)
1326  BL flux of atmospheric tracer 25 (conventionally SOOT)
1327  BL flux of atmospheric tracer 26 (conventionally VOLCANIC ASH)
1328  BL flux of atmospheric tracer 27 (conventionally SULPHURIC ACID)
1329  BL flux of atmospheric tracer 28 (conventionally AMMONIUM SULPHAT
1330  BL flux of atmospheric tracer 29 (conventionally MINERAL)
1331  BL flux of atmospheric total aerosol
1332  Ice Possible (0- -20 Celcius, RH =>70%) 1.0=possible;0.0 not
1333  Liquid water in any region averaged during all-sky conditions
1334  q*w  (specific humidity * vertical velocity)
1335  Fraction of time pressure level is above model surface
1336  U * geopotential height
1337  V * geopotential height
1338  Virtual temperature
1339  Saturation mixing ratio
1340  Reserved for VAR project : B.Ingleby/M Wlasek
1369  Reserved for VAR project : B Ingleby/M Wlasek
1370  SO4 aerosol: Aitken mode. Units: mass mixing ratio (kg/kg)
1371  SO4 aerosol: accumulation mode. Units: mass mixing ratio (kg/kg)
1372  SO4 aerosol: dissolved mode. Units: mass mixing ratio (kg/kg)
1373  SO4 aerosol: DMS mass mixing ratio
1374  SO4 aerosol: SO2 mass mixing ratio
1375  SO4 Methyl Sulphonic Acid mass mixing ratio kg/kg
1376  SO4 aerosol: column mass loading (kg/m2)
1377  SO4 aerosol: SW radiative forcing (W/m2)
1379  NH3 gas mass mass mix ratio kg/kg
1380  Van-Genuchten 'B' parameter
1381  Clapp-Hornberger 'B' exponent
1382  Leaf area index of vegetated fraction
1383  Canopy height of vegetated fraction
1384  Canopy conductance
1385  Unfrozen soil moisture fraction
1386  Frozen soil moisture fraction
1387  Transpiration
1388  Gross Primary Productivity
1389  Net Primary Productivity
1390  Plant Respiration
1391  Fractional covering of functional types
1392  LAI of veg functional types
1393  Canopy height of vegetation functional types
1394  Disturbed fraction of vegetation
1395  Soil albedo
1396  Snow soot
1397  Soil carbon content
1398  Accumulated net primary productivity on tiles (kg C/m2/sec)
1399  Mountain torque per unit area
1400  Surface dry deposition flux of SO2 kg/m2/s
1401  Surface dry deposition flux of SO4 Aitken mode kg/m2/s
1402  Surface dry deposition flux of SO4 ACCUMULATION mode kg/m2/s
1403  Surface dry deposition flux of SO4 DISSOLVED mode kg/m2/s
1404  RESIST_B for SO2 (Note: 1400-1416 are Sulphur Cycle BL diagnostic
1405  RESIST_B for SO2 Aitken mode
1406  RESIST_B for SO2 ACCUMULATION  mode
1407  RESIST_B for SO2 DISSOLVED mode kg/m2/s
1408  RESIST_S for SO2
1409  RESIST_S for SO2 Aitken mode
1410  RESIST_S for SO2 ACCUMULATION mode
1411  RESIST_S for SO2 DISSOLVED  mode
1412  Dry deposition velocity for SO2
1413  Dry deposition velocity for SO4 Aitken mode
1414  Dry deposition velocity for SO4 ACCUMULATION mode
1415  Dry deposition velocity for SO4 DISSOLVED mode
1416  Aerodynamic resistance 1/CDSTD after TSTEP
1417  SO2 scavenged by convective precipitation kg(S)/m2/tstep
1418  SO4 Aitken scavenged by convective precipitation kg(S)/m2/tstep
1419  SO4 accumulation scavenged by convective precipitation kg(S)/m2/t
1420  SO4 dissolved scavenged by convective precipitation kg(S)/m2/tste
1421  SO2 scavenged by large scale precipitation kg(S)/m2/tstep
1422  SO4 Aitken scavenged by large scale precipitation kg(S)/m2/tstep
1423  SO4 accumulation scavenged by large scale precipitation kg(S)/m2/
1424  SO4 dissolved scavenged by large scale precipitation kg(S)/m2/tst
1425  Layer liquid cloud amount in layers
1426  Layer frozen cloud amount in layers
1427  Layer cloud frequency in each layer
1428  Net energy change this period J/m**2
1429  Thickness (alternative code for other packages eg without BRLEV s
1430  Thickness tendency (as above)
1431  du/dp
1432  dv/dp
1433  dtheta/dp
1434  dtheta_e/dp
1435  Moist PV
1436  Magnitude Geostrophic Deformation
1437  EW cpt, Geostrophic deformation axis
1438  NS cpt, Geostrophic deformation axis
1439  Magnitude Grad(ThetaW)
1440  EW cpt Grad(ThetaW)
1441  NS cpt Grad(ThetaW)
1442  Geostrophic Deformation wrt ThetaW
1443  Geostrophic Relative Vorticity
1444  Magnitude Grad(Theta)
1445  EW cpt Grad(Theta)
1446  NS cpt Grad(Theta)
1447  Omega eqn. inversion with diabatic forcing
1448  Omega eqn. inversion with divQ & diabatic forcing
1449  Omega eqn. inversion with divQ only
1450  Omega eqn. inversion with divQ 800-600mb
1451  Omega eqn. inversion with divQ 600-100mb
1452  Omega eqn. inversion with divQ 1000-800mb
1453  DelSqd(ThetaW)
1454  Local value of magnitude of grad(ThetaW)
1455  Theta Frontal Speed - advection of magnitude (theta_w)
1456  Advection of theta_w
1457  Geostrophic advection of theta_w
1458  DelSqd(Theta)
1459  Local value of magnitude of grad(Theta)
1460  Theta Frontal Speed - advection of magnitude (theta)
1461  Advection of theta
1462  Geostrophic advection of theta
1463  Ambient noise
1464  Droplet number conc * cloud amount
1465  Layer cloud lwc * cloud amount
1466  SO4 ccn mass conc * cond samp weight
1467  Conditional sampling weight
1468  2-D effective radius * 2-D re weight
1469  Weight for 2-D effective radius
1470  Advection of mixing ratio
1471  Advection of dew point temperature
1472  Advection of relative humidity
1473  Advection of relative vorticity
1474  Advection of absolute vorticity
1475  'X' Component of Isallobaric Wind (corrected from u July 1998)
1476  'Y' Component of Isallobaric Wind (corrected from v July 1998)
1477  Scalar Divergence
1478  Showalter index (an instability index)
1479  Total total index (an instability index)
1480  Sweat index (an instability index)
1481  Lifted index (an instability index)
1482  Cape (convective available potential energy)
1483  k index (an instability index)
1484  Saturation vapour pressure
1485  Modified refractivity index, M (M units)
1486  Vertical gradient of modified refractivity, dM/dz (M/km)
1487  Minimum of dM/dz with height (M/km)
1488  Maximum of dM/dz with height (M/km)
1489  Radar duct intensity (max wavelength trapped) (cm)
1490  Refractivity index, N (N units)
1491  Fresh soot mass mix ratio kg/kg
1492  Aged soot mass mix ratio kg/kg
1493  Soot in cloud mix ratio kg/kg
1494  SO4 Aitken surface settlement flux kg/m2/s
1495  SO4 Accumulation surface settlement flux kg/m2/s
1496  MSA (methyl sulphonic acid) mass mixing ratio flux kg/kg/s
1497  NH3 (ammonia) depletion after tstep kg/kg
1498  Ship noise (from passing ships)
1499  Leaf turnover rate of plant functional types (/360days)
1500  Accumulated leaf turnover rate on PFTs (1/sec)
1501  Accumulated phenological leaf turnover rate on PFTs (sec-1)
1502  Accumulated wood respiration on PFTs (kg C/m2/sec)
1503  Accumulated soil respiration (kg C/m2/sec)
1504  Canopy & surface water content on tiles (kg/m2)
1505  Canopy & surface capacity on tiles (kg/m2)
1506  Infiltration enhancement factor on tiles
1507  Snow grain size (microns)
1508  Snow temperature (K)
1509  Radiative surface temperature (K)
1510  Surface temperature on tiles (K)
1511  Roughness length on tiles (m)
1512  Vegetation carbon on Plant Functional Types (PFTs) (kg C / m2)
1513  Gridbox mean vegetation carbon  (kg C / m2)
1514  Phenological leaf turnover rate on PFTs (1/sec)
1515  Litter carbon on Plant Functional Types (PFTs) (kg C / m2 / year)
1516  Gridbox mean litter carbon (kg C / m2 / year)
1517  Canopy evaporation on tiles (kg / m2 /sec)
1518  Evapotranspiration from the soil (kg / m2 /sec)
1519  Gross primary productivity on plant functional types (kg C/m2/sec
1520  Sensible heat flux on tiles (W/m2)
1521  Net primary productivity on plant functional types (kg C/m2/sec)
1522  Plant respiration on plant functional types (kg C/m2/sec)
1523  Soil respiration (kg C/m2/sec)
1524  Bulk Richardson number on tiles
1525  Fractional snow cover
1526  Rate of evaporation from soil surface  kg/m2/s
1527  Rate of evaporation from canopy  kg/m2/s
1528  Rate of sublimation from surface (gridbox mean) kg/m2/s
1529  Rate of transpiration  kg/m2/s
1530  Rate of snow melt on land kg/m2/s
1531  Rate of canopy throughfall kg/m2/s
1532  Rate of surface runoff kg/m2/s
1533  Rate of sub-surface runoff kg/m2/s
1534  Turbulent mixing ht after Boundary Layer m
1535  Stable BL indicator
1536  Stratocumulus over stable BL indicator
1537  Well_mixed BL indicator
1538  Decoupled SC not over CU indicator
1539  Decoupled SC over CU indicator
1540  Cumulus-capped BL indicator
1541  NH3 Surface dry deposition flux kg/m2/s
1542  Fresh soot surface dry deposition flux kg/m2/s
1543  Aged soot surface dry deposition flux kg/m2/s
1544  Soot in cloud occult deposition flux kg/m2/s
1545  NH3 scavenged by large scale precipitation kg/m2/s
1546  Soot scavenged by large scale precipitation kg/m2/s
1547  SO2 scavenged by large scale precipitation kg/m2/s
1548  SO4 Aitken scavenged by large scale precipitation kg/m2/s
1549  SO4 accumulation scavenged by large scale precipitation kg/m2/s
1550  SO4 dissolved scavenged by large scale precipitation kg/m2/s
1551  NH3 scavenged by convective precipitation kg/m2/s
1552  Soot scavenged by convective precipitation kg/m2/s
1553  SO2 scavenged by convective precipitation kg/m2/s
1554  SO4 Aitken scavenged by convective precipitation kg/m2/s
1555  SO4 accumulation scavenged by convective precipitation kg/m2/s
1556  SO4 dissolved scavenged by convective precipitation kg/m2/s
1557  Potential evaporation amount kg/m2/tstep
1558  Potential evaporation rate kg/m2/s
1559  Soil moisture availability factor
1560  CO2 ocean flux kg/m**2/s
1561  CO2 surface emissions kg/m**2/s
1562  CO2 land surface flux kg/m**2/s
1563  CO2 total flux to atmosphere kg/m**2/s
1564  Precipitation rate codes for symbol plotting
1565  Present Weather codes (ww-code)
1566  T before dynamics
1567  T after dynamics
1568  q before dynamics
1569  q after dynamics
1570  x component of E vector = (v'\*\*2 - u'\*\*2)
1571  y component of E vector = -u'v'
1572  z component of E vector - proportional to v'T'
1573  Supercooled liquid water content
1574  Supercooled rain
1575  Rain fraction
1576  Cloud composite
1577  Accumulated precipitation for 6 hours (for fieldsfiles, 1998)
1578  Modified maximum gravity wave stress magnitude. Units N/m**2
1579  U in the GWD surface layer
1580  V in the GWD surface layer
1581  N in the GWD surface layer
1582  GWD surface Froude number
1583  GWD Blocked Layer Depth
1584  Percent of hydro GWD that is linear
1585  Percent of time with Hydr. Jumps
1586  Percent of time with Lee Waves
1587  Percent of time with blocked flow
1588  x component of GW saturation stress
1589  y component of GW saturation stress
1590  x component of GW jump stress
1591  y component of GW jump stress
1592  x component of GW wake stress
1593  y component of GW wake stress
1594  x component of GW lee stress
1595  y component of GW lee stress
1596  u-accel from GWD blocked flow
1597  v-accel from GWD blocked flow
1598  Shallow convection indicator
1599  Cumulus over orography indicator
1600  Net surface water flux/kg m**-2 s**-1
1601  Pressure tendency (for Thailand, 1997, system has no LBPROC)
1602  'X' Component of Thermal Wind (corrected from u July 1998)
1603  'Y' Component of Thermal Wind (corrected from v July 1998)
1604  Saturated specific humidity
1605  Isallobaric direction ( true )
1606  Isallobaric wind speed
1607  Geostrophic wind direction ( true )
1608  Geostrophic wind speed
1609  Thermal wind direction ( true )
1610  ISCCP C2 high thin cloud amount
1611  ISCCP C2 high medium cloud amount
1612  ISCCP C2 high thick cloud amount
1613  ISCCP C2 mid-level thin cloud amount
1614  ISCCP C2 mid-level thick cloud amount
1615  ISCCP C2 low thin cloud amount
1616  ISCCP C2 low thick cloud amount
1617  Sea Surface Temperature (SST) minus Dew point (for forecasting se
1618  500-850 Potential Instability Index (measure of thunder risk)
1619  Boyden Index
1620  Height of 2 deg C isotherm in feet (Navy requirement)
1621  Casswell maximum vertical velocity (for NMC)
1622  Casswell height of maximum vertical velocity (for NMC)
1623  Pseudo WV image from Horace - estimated radiance received
1624  Pseudo IR image from Horace - estimated radiance received
1625  Frequency unit type: 1=Hz (eg for noise fields)
1626  Areal fraction of intensive grassland
1627  Areal fraction of extensive grassland
1628  SHEAR-DOMINATED UNSTABLE BL INDICATOR
1629  CIN Convective INhibition (units J/kg)
1630  Soil clay fraction (ancil)
1631  Soil silt fraction (ancil)
1632  Soil sand fraction (ancil)
1633  Soil mass fraction, dust divisions 1->6
1634  Mineral dust mixing ratio divisions 1->6
1635  Dust emission flux (kg m-2 s-1), divisions 1->6
1636  Threshold friction vel. dust divisions 1->6
1637  Dry threshold friction vel. dust divisions 1->6
1638  Friction velocity for dust
1639  Dust dry dep flux from 1st lev, divisions 1->6
1640  Dust dry dep flux from 2nd lev, divisions 1->6
1641  Dust wet dep flux from LSP, divisions 1->6
1642  Dust wet dep flux from CONV, divisions 1->6
1643  Surface dust concentration.
1644  Dust concentration between 2000-5000ft.
1645  TOTAL DUST CONC (microg/m3)
1646  Dust emission fraction
1647  Mineral Dust Optical Code.
1648  RESERVED FOR DAVID WALTERS 8.1
1649  VIS AT 1.5M  (incl ppn and dust)
1650  VIS IN DUST ONLY AT 1.5M
1651  MW surface emissivity
1652  MW surface emissivity (stdev)
1653  Reserved Mineral Dust numbers for S Woodward
1654  Reserved Mineral Dust numbers for S Woodward
1655  Reserved Mineral Dust numbers for S Woodward
1656  Reserved Mineral Dust numbers for S Woodward
1657  Reserved Mineral Dust numbers for S Woodward
1658  Reserved Mineral Dust numbers for S Woodward
1659  Reserved Mineral Dust numbers for S Woodward
1660  Reserved Mineral Dust numbers for S Woodward
1661  Reserved Mineral Dust numbers for S Woodward
1662  Reserved Mineral Dust numbers for S Woodward
1663  Reserved Mineral Dust numbers for S Woodward
1664  Reserved Mineral Dust numbers for S Woodward
1665  Reserved Mineral Dust numbers for S Woodward
1666  Reserved Mineral Dust numbers for S Woodward
1667  Reserved Mineral Dust numbers for S Woodward
1668  Reserved Mineral Dust numbers for S Woodward
1669  Reserved Mineral Dust numbers for S Woodward
1670  Reserved Mineral Dust numbers for S Woodward
1671  Reserved Mineral Dust numbers for S Woodward
1672  Reserved Mineral Dust numbers for S Woodward
1673  Reserved Mineral Dust numbers for S Woodward
1674  Reserved Mineral Dust numbers for S Woodward
1675  Reserved Mineral Dust numbers for S Woodward
1676  Reserved Mineral Dust numbers for S Woodward
1677  Reserved Mineral Dust numbers for S Woodward
1678  Reserved Mineral Dust numbers for S Woodward
1679  Reserved Mineral Dust numbers for S Woodward
1680  Reserved Mineral Dust numbers for S Woodward
1681  Rain water [kg/kg]
1682  Cape inhibition
1683  Mass mixing ratio fresh smoke
1684  Mass mixing ratio aged smoke
1685  Mass mixing ratio cloud smoke
1686  Dry deposition fresh smoke
1687  Dry deposition aged smoke
1688  Occult deposition cloud smoke
1689  Frozen cloud water (pristine crystals) [kg/kg]
1690  C.A.T. probability derived from fieldsfiles, not model
1691  LS rainout of biomass smoke
1692  LS washout of biomass smoke
1693  Convective scavenging of smoke
1694  Maximum wind gust at 10m
1695  Obukhov length
1696  Frictional velocity
1697  Nominal 3D convective rainfall rate (kg/m2/s)
1698  Nominal 3D convective snowfall rate (kg/m2/s)
1699  Boundary Layer Thermal Strength (m/s)
1701  OD
1702  OP           'STOCHEM' TROP. CHEMISTRY MODEL:
1703  OH       A fuller list of descriptions is expected from
1704  NO          Colin Johnson
1705  NO2      12/8/99
1706  NO3
1707  N2O5
1708  CO
1709  CH4
1710  HCHO
1711  O3
1712  H2
1713  HNO3
1714  H2O2
1715  CH3O2
1716  HO2
1717  C2H6
1718  C2H5O2
1719  CH3CHO
1720  CH3COO2
1721  PAN
1722  CH3OOH
1723  NC4H10
1724  SC4H9O2
1725  CH3COE
1726  C2H4
1727  C3H6
1728  C3H8
1729  C3H7O2
1730  C3H7OOH
1731  C2H5OOH
1732  C4H9OOH
1733  CH3OH
1734  ACETONE
1735  ACETO2
1736  CH3COX
1737  CH2O2C
1738  MGLYOX
1739  CH3CHX
1740  C5H8
1741  RO2IP1
1742  MVK
1743  RO2IP2
1744  ISOPOOH
1745  MVKOOH
1746  RNC2H4
1747  RNC3H6
1748  RNC5H8
1749  NAER
1750  HO2NO2
1751  H2O
1752  EXTRA
1801  CCA FROM DEEP CONVECTION (no dimension)
1802  CCA FROM MID-LEVEL CONVECTION (no dimension)
1803  CCA FROM SHALLOW CONVECTION (no dimension)
1860  Volumetric Mixing Ratio
1861  Mass Mixing Ratio
1862  Particle Number Density [Particles per molecule of air]
1863  Particle Number Density [Number per cm^3]
1864  Particle Area Density [m^2 per m^3]
1865  Molar Density [mol per m^3]
1866  Molecular Density [molecules per cm^3]
1867  Mass Density [kg per m^3]
1868  Molarity [mol per litre]
1869  Molality [mol per kg]
1870  Molecular flux density [molecules per (s.cm^3)]
1871  Molar flux density [mol per (s.m^3)]
1872  Mass flux density [kg per (s.m^3)]
1873  Particle flux density [Number per (s.m^3)]
1874  Molecular surface flux [molecules per (s.cm^2)]
1875  Molar surface flux [mol per (s.m^2)]
1876  Mass surface flux [kg per (s.m^2)]
1877  Particle surface flux [Number per (s.m^2)]
1878  Mass flux on atmosphere theta levels (kg m-2 s-1)
1879  THERMAL ROUGHNESS LENGTH ON TILES          M
1881  WIND-SEA ENERGY AT PEAK FREQ  (m^2/Hz/rad)
1882  SWELL ENERGY AT PEAK FREQ    (m^2/Hz/rad)
1883  2ND SWELL ENERGY AT PEAK FREQ  (m^2/Hz/rad)
1884  3RD SWELL ENERGY AT PEAK FREQ   (m^2/Hz/rad)
1885  WIND-SEA MEAN DIRECTIONAL SPREAD  (deg)
1886  SWELL MEAN DIRECTIONAL SPREAD    (deg)
1887  2ND SWELL MEAN DIRECTIONAL SPREAD(deg)
1888  3RD SWELL MEAN DIRECTIONAL SPREAD (deg)
1890  Additional radiation information between timesteps
1891  Total downward PAR flux at the surface
1892  Direct component of PAR flux at surface
1893  Surface Tile Fractions
1894  Leaf Area Indices on PFTS
1895  Canopy Heights on PFTS
1896  Stomatal Conductance on PFTS (m/s)
1900  Magnitude of wind stress
1901  River outflow (Kg/m2/s)
1902  River water storage (Kg)
1903  Total gridbox river inflow (Kg/s)
1904  Total gridbox river outflow (Kg/s)
1905  Reserved River Routing for Cyndy Bunton.
1906  INLAND BASIN OUTFLOW ON TRIP RIVER ROUTING GRID    KG/M2/S
1907  INLAND BASIN OUTFLOW ON ATMOSPHERE GRID    KG/M2/S
1908  Reserved River Routing for Cyndy Bunton.
1909  Reserved River Routing for Cyndy Bunton.
1910  Reserved River Routing for Cyndy Bunton.
1911  Reserved River Routing for Cyndy Bunton.
1912  Reserved River Routing for Cyndy Bunton.
1913  Reserved River Routing for Cyndy Bunton.
1914  Reserved River Routing for Cyndy Bunton.
1915  Reserved River Routing for Cyndy Bunton.
1916  Reserved River Routing for Cyndy Bunton.
1917  Reserved River Routing for Cyndy Bunton.
1918  Reserved River Routing for Cyndy Bunton.
1919  Reserved River Routing for Cyndy Bunton.
1920  Reserved River Routing for Cyndy Bunton.
1921  Graupel (kg/kg)
1922  General Atmos Codes
1923  General Atmos Codes
1924  General Atmos Codes
1925  General Atmos Codes
1926  General Atmos Codes
1927  General Atmos Codes
1928  General Atmos Codes
1929  TOTAL SPECIFIC HUMIDITY DA INCREMENT
1930  Ozone Tropopause Level
1931  Ozone Tropopause Height
1932  Thermal Tropopause Level
1933  Thermal Tropopause Height
1934  General Atmos Codes
1935  General Atmos Codes
1936  General Atmos Codes
1937  General Atmos Codes
1938  General Atmos Codes
1939  General Atmos Codes
1940  General Atmos Codes
1941  General Atmos Codes
1942  General Atmos Codes
1943  General Atmos Codes
1944  General Atmos Codes
1945  General Atmos Codes
1946  General Atmos Codes
1947  General Atmos Codes
1948  General Atmos Codes
1949  General Atmos Codes
1950  Lumb Snow
1951  Topmelt for individual ice categories (W m-2)
1952  Botmelt for individual ice categories  (W m-2)
1953  Fractional ice cover for individual ice categories
1954  Grid box mean ice depth for individual ice categories (m)
1955  Local snow depth over ice for individual ice categories (m)
1956  Rate of change of sea-ice concentration for individual ice categories (s-1)
1957  Rate of change of sea-ice depth for individual ice categories (m s-1)
1958  Rate of change of snowdepth on sea-ice for individual ice categories (m s-1)
1959  Grid box mean snowdepth for individual ice categories (m)
1960  Wave/ Swell peak period
1961  Wave/ Wavetrain wave height
1962  Wave/ Wavetrain upcrossing period
1963  Wave/ Wavetrain mean direction
1964  Wave/ Mean energy at a given frequency
1965  Wave/ Mean direction at a given frequency
1966  Wave/ Total mean direction
1967  Wave/ Total mean period
1968  Wave/ Windsea mean directon
1969  Wave/ Windsea mean period
1970  Wave/ Swell mean direction
1971  Wave/ Swell mean period
1972  Wave/ Energy in each frequency (1d spectrum)
1973  Wave/ Principal direction of 1d spectrum
1974  Wave/ Mean direction of 1d spectrum
1975  Wave/ 2d energy spectrum
1976  reserved for any further Wave
1977  reserved for any further Wave
1978  reserved for any further Wave
1979  reserved for any further Wave
1980  reserved for any further Wave
1981  WAFC Mean Icing Potential (Index)
1982  WAFC Max. Icing Potential (Index)
1983  WAFC Mean In-Cloud Turbulence Potential (Index)
1984  WAFC Max. In-Cloud Turbulence Potential (Index)
1985  WAFC Mean CAT Potential (Index)
1986  WAFC Max. CAT Potential (Index)
1987  WAFC CB Horizontal Extent (Index)
1988  WAFC Pressure at CB Base (Pa)
1989  WAFC Pressure at CB Top (Pa)
1990  WAFC Pressure at Embedded CB Base (Pa)
1991  WAFC Pressure at Embedded CB Top (Pa)
1992  WAFC ICAO Height at CB Base (kft)
1993  WAFC ICAO Height at CB Top (kft)
1994  WAFC ICAO Height at Embedded CB Base (kft)
1995  WAFC ICAO Height at Embedded CB Top (kft)
1996  Cloud Top Height at Mean Icing Potential (m)
2001  -70 deg.C level
2002  Reserved
2003  Reserved
2004  Reserved
2005  Reserved
2006  Reserved
2007  Reserved
2008  Reserved
2009  Reserved
2010  Reserved
2011  Soil Bulk Density (kg/m^3 ).. reserved for D Walters.
2020  PM10 mass concentration (microgrammes per cubic metre)
2021  PM2.5 mass concentration (microgrammes per cubic metre)
2022  Sulphate contribution to PM10 concentration (microgrammes per cubic metre)
2023  Sulphate contribution to PM2.5 concentration (microgrammes per cubic metre)
2024  Black carbon contribution to PM10 concentration (microgrammes per cubic metre)
2025  Black carbon contribution to PM2.5 concentration (microgrammes per cubic metre)
2026  Biomass burning contribution to PM10 concentration (microgrammes per cubic metre)
2027  Biomass burning contribution to PM2.5 concentration (microgrammes per cubic metre)
2028  OCFF contribution to PM10 concentration (microgrammes per cubic metre)
2029  OCFF contribution to PM2.5 concentration (microgrammes per cubic metre)
2030  SOA contribution to PM10 concentration (microgrammes per cubic metre)
2031  SOA contribution to PM2.5 concentration (microgrammes per cubic metre)
2032  Sea salt contribution to PM10 concentration (microgrammes per cubic metre)
2033  Sea salt contribution to PM2.5 concentration (microgrammes per cubic metre)
2034  Dust contribution to PM10 concentration (microgrammes per cubic metre)
2035  Dust contribution to PM2.5 concentration (microgrammes per cubic metre)
2036  Nitrate contribution to PM10 concentration (microgrammes per cubic metre)
2037  Nitrate contribution to PM2.5 concentration (microgrammes per cubic metre)
2038  TURBULENT KINETIC ENERGY, m**2 s**-2
2039  SELF COVARIANCE OF THETAL', K**2
2040  SELF COVARIANCE OF QW', kg**2 kg**-2
2041  CORRELATION OF THETAL' AND QW', K kg kg**-1
2042  COUNTER GRADIENT TERM OF TAUX, kg m**-1 s**-2
2043  COUNTER GRADIENT TERM OF TAUY, kg m**-1 s**-2
2044  COUNTER GRADIENT TERM OF FTL, W m**-2
2045  COUNTER GRADIENT TERM OF FQW, W m**-2
2046  MIXING LENGTH, m
2047  PRODUCTION RATE OF TKE BY SHEAR, m**2 s**-3
2048  PRODUCTION RATE OF TKE BY BUOYANCY, m**2 s**-3
2049  DISSIPATION RATE OF TKE, m**2 s**-3
2050  NON-DIM DIFFUSION COEFS FOR UV (SM) , (no dimension)
2051  NON-DIM DIFFUSION COEFS FOR TQ (SH) , (no dimension)
2052  NON-GRADIENT BUOYANCY FLUX, m**2 s**-3
2053  CLOUD FRACTION IN THE TKE SCHEMES, (no dimension)
2054  CONDENSED WATER IN THE TKE SCHEMES, kg kg**-1
2055  STD. DEV. OF GAUSSIAN, kg kg**-1
2056  Height of mixed layer to evaluate the non-gradient buoyancy flux (m)
2057  Mean shear Clear Air Turbulence (% probability of occurrence)
2058  Max shear Clear Air Turbulence  (% probability of occurrence)
2059  Mean Mountain/gravity wave stress on certain pressure levels (n/m**2)
2060  Max Mountain/gravity wave stress on certain pressure levels  (n/m**2)
2061  SFERICS (lightning 'count' per area per time) (See Note d above).
2062  Helicopter-Triggered Lightning Risk Index (0 - 30). (See Note e above).
2063  Lightning Index (0, 1, 10)'. (See Note f above).
2064  Lightning Risk (within 50km) Index (5 - 1)'.(See Note f above).
2065  CONVECTIVE RAIN RATE FOR NAME
2066  CONVECTIVE SNOW RATE FOR NAME
2101  O3 MASS MIXING RATIO                 (kg/kg)
2102  NO MASS MIXING RATIO                 (kg/kg)
2103  NO3 MASS MIXING RATIO                (kg/kg)
2104  NO2 MASS MIXING RATIO                (kg/kg)
2105  N2O5 MASS MIXING RATIO               (kg/kg)
2106  HO2NO2 MASS MIXING RATIO             (kg/kg)
2107  HONO2 MASS MIXING RATIO              (kg/kg)
2108  H2O2 MASS MIXING RATIO               (kg/kg)
2109  CH4 MASS MIXING RATIO                (kg/kg)
2110  CO MASS MIXING RATIO                 (kg/kg)
2111  HCHO MASS MIXING RATIO               (kg/kg)
2112  MeOOH MASS MIXING RATIO              (kg/kg)
2113  HONO MASS MIXING RATIO               (kg/kg)
2114  C2H6 MASS MIXING RATIO               (kg/kg)
2115  EtOOH MASS MIXING RATIO              (kg/kg)
2116  MeCHO MASS MIXING RATIO              (kg/kg)
2117  PAN MASS MIXING RATIO                (kg/kg)
2118  C3H8 MASS MIXING RATIO               (kg/kg)
2119  n-PrOOH MASS MIXING RATIO            (kg/kg)
2120  i-PrOOH MASS MIXING RATIO            (kg/kg)
2121  EtCHO MASS MIXING RATIO              (kg/kg)
2122  Me2CO MASS MIXING RATIO              (kg/kg)
2123  MeCOCH2OOH MASS MIXING RATIO         (kg/kg)
2124  PPAN MASS MIXING RATIO               (kg/kg)
2125  MeONO2 MASS MIXING RATIO             (kg/kg)
2126  O3S MASS MIXING RATIO                (kg/kg)
2127  C5H8 MASS MIXING RATIO               (kg/kg)
2128  ISOOH MASS MIXING RATIO              (kg/kg)
2129  ISON MASS MIXING RATIO               (kg/kg)
2130  MACR MASS MIXING RATIO               (kg/kg)
2131  MACROOH MASS MIXING RATIO            (kg/kg)
2132  MPAN MASS MIXING RATIO               (kg/kg)
2133  HACET MASS MIXING RATIO              (kg/kg)
2134  MGLY MASS MIXING RATIO               (kg/kg)
2135  NALD MASS MIXING RATIO               (kg/kg)
2136  HCOOH MASS MIXING RATIO              (kg/kg)
2137  MeCO3H MASS MIXING RATIO             (kg/kg)
2138  MeCO2H MASS MIXING RATIO             (kg/kg)
2140  ISO2 MASS MIXING RATIO               (kg/kg)
2141  Cl MASS MIXING RATIO                 (kg/kg)
2142  ClO MASS MIXING RATIO                (kg/kg)
2143  Cl2O2 MASS MIXING RATIO              (kg/kg)
2144  OClO MASS MIXING RATIO               (kg/kg)
2145  Br MASS MIXING RATIO                 (kg/kg)
2147  BrCl  MASS MIXING RATIO              (kg/kg)
2148  BrONO2 MASS MIXING RATIO             (kg/kg)
2149  N2O MASS MIXING RATIO                (kg/kg)
2151  HOCl MASS MIXING RATIO               (kg/kg)
2152  HBr MASS MIXING RATIO                (kg/kg)
2153  HOBr MASS MIXING RATIO               (kg/kg)
2154  ClONO2 MASS MIXING RATIO             (kg/kg)
2155  CFCl3 MASS MIXING RATIO              (kg/kg)
2156  CF2Cl2 MASS MIXING RATIO             (kg/kg)
2157  MeBr MASS MIXING RATIO               (kg/kg)
2158  N MASS MIXING RATIO                  (kg/kg)
2159  O3P MASS MIXING RATIO                (kg/kg)
2160  MACRO2 MASS MIXING RATIO             (kg/kg)
2170  H2 MASS MIXING RATIO                 (kg/kg)
2171  DMS MASS MIXING RATIO                (kg/kg)
2172  SO2 MASS MIXING RATIO                (kg/kg)
2173  H2SO4 MASS MIXING RATIO              (kg/kg)
2174  MSA MASS MIXING RATIO                (kg/kg)
2175  DMSO MASS MIXING RATIO               (kg/kg)
2176  NH3 MASS MIXING RATIO                (kg/kg)
2177  CS2 MASS MIXING RATIO                (kg/kg)
2178  COS MASS MIXING RATIO                (kg/kg)
2179  H2S MASS MIXING RATIO                (kg/kg)
2180  H MASS MIXING RATIO                  (kg/kg)
2181  OH MASS MIXING RATIO                 (kg/kg)
2182  HO2 MASS MIXING RATIO                (kg/kg)
2183  MeOO MASS MIXING RATIO               (kg/kg)
2184  EtOO MASS MIXING RATIO               (kg/kg)
2185  MeCO3 MASS MIXING RATIO              (kg/kg)
2186  n-PrOO MASS MIXING RATIO             (kg/kg)
2187  i-PrOO MASS MIXING RATIO             (kg/kg)
2188  EtCO3 MASS MIXING RATIO              (kg/kg)
2189  MeCOCH2OO MASS MIXING RATIO          (kg/kg)
2190  MeOH MASS MIXING RATIO               (kg/kg)
2191  MONOTERPENE MASS MIXING RATIO        (kg/kg)
2192  SEC_ORG MASS MIXING RATIO            (kg/kg)
2193  C3H6 MASS MIXING RATIO               (kg/kg)
2194  SO3 MASS MIXING RATIO                (kg/kg)
2195  C4H9OOH MASS MIXING RATIO            (kg/kg)
2196  MEK MASS MIXING RATIO                (kg/kg)
2197  TOLUENE MASS MIXING RATIO            (kg/kg)
2198  LUMPED N (as NO2) MMR                (kg/kg)
2199  LUMPED Br (as BrO) MMR               (kg/kg)
2200  LUMPED Cl (as HCl) MMR               (kg/kg)
2201  NUCLEATION MODE (SOLUBLE) NUMBER     (Particles per molecule of air)
2202  NUCLEATION MODE (SOLUBLE) H2SO4 MMR  (kg/kg)
2203  AITKEN MODE (SOLUBLE) NUMBER         (Particles per molecule of air)
2204  AITKEN MODE (SOLUBLE) H2SO4 MMR      (kg/kg)
2205  AITKEN MODE (SOLUBLE) BC MMR         (kg/kg)
2206  AITKEN MODE (SOLUBLE) OC MMR         (kg/kg)
2207  ACCUMULATION MODE (SOLUBLE) NUMBER   (Particles per molecule of air)
2208  ACCUMULATION MODE (SOL) H2SO4 MMR    (kg/kg)
2209  ACCUMULATION MODE (SOL) BC MMR       (kg/kg)
2210  ACCUMULATION MODE (SOL) OC MMR       (kg/kg)
2211  ACCUMULATION MODE (SOL) SEA SALT MMR (kg/kg)
2212  ACCUMULATION MODE (SOL) DUST MMR     (kg/kg)
2213  COARSE MODE (SOLUBLE) NUMBER         (Particles per molecule of air)
2214  COARSE MODE (SOLUBLE) H2SO4 MMR      (kg/kg)
2215  COARSE MODE (SOLUBLE) BC MMR         (kg/kg)
2216  COARSE MODE (SOLUBLE) OC MMR         (kg/kg)
2217  COARSE MODE (SOLUBLE) SEA SALT MMR   (kg/kg)
2218  COARSE MODE (SOLUBLE) DUST MMR       (kg/kg)
2219  AITKEN MODE (INSOLUBLE) NUMBER       (Particles per molecule of air)
2220  AITKEN MODE (INSOLUBLE) BC MMR       (kg/kg)
2221  AITKEN MODE (INSOLUBLE) OC MMR       (kg/kg)
2222  ACCUMULATION MODE (INSOLUBLE) NUMBER (Particles per molecule of air)
2223  ACCUMULATION MODE (INSOLUBLE) DUST   (kg/kg)
2224  COARSE MODE (INSOLUBLE) NUMBER       (Particles per molecule of air)
2225  COARSE MODE (INSOLUBLE) DUST MMR     (kg/kg)
2226  NUCLEATION MODE (SOLUBLE) OC  MMR    (kg/kg)
2227  AITKEN MODE (SOLUBLE) SEA SALT MMR   (kg/kg)
2228  NUCLEATION MODE (SOLUBLE) OC2 MMR    (kg/kg)
2229  AITKEN MODE (SOLUBLE) OC2 MMR        (kg/kg)
2230  ACCUMULATION MODE (SOLUBLE) OC2 MMR  (kg/kg)
2231  COARSE MODE (SOLUBLE) OC2 MMR        (kg/kg)
2232  NUCLEATION MODE (SOLUBLE) NH4 MMR    (kg/kg)
2233  AITKEN MODE (SOLUBLE) NH4 MMR        (kg/kg)
2234  ACCUMULATION MODE (SOLUBLE) NH4 MMR  (kg/kg)
2235  COARSE MODE (SOLUBLE) NH4 MMR        (kg/kg)
2236  NUCLEATION MODE (SOLUBLE) NO3 MMR    (kg/kg)
2237  AITKEN MODE (SOLUBLE) NO3 MMR        (kg/kg)
2238  ACCUMULATION MODE (SOLUBLE) NO3 MMR  (kg/kg)
2239  COARSE MODE (SOLUBLE) NO3 MMR        (kg/kg)
2249  PASSIVE O3 MASS MIXING RATIO         (kg/kg)
2250  AGE OF AIR IN SECONDS                (s)
3845  Photolysis rate JO2 (s-1)
3846  Photolysis rate JO3P (s-1)
====  ==================================================================================

Notes
=====

   a. Field codes 30 to 34 are used to indicate cloud AMOUNTS by
      coding zero for the vertical co-ordinate type.  When a
      vertical co-ordinate type in the range 135 to 138 is coded,
      these field codes represent the upper or lower boundaries of
      clouds or contrails.

   b. Some of the original M20 fields were duplicates of existing fields and
      have been removed, others have been re-defined to conform to official PP
      standards.

   c. Type I clear-sky quantities are those from those points
      and times when no clouds are in fact present.
      Type II are those that would be produced everywhere
      in the absence of clouds but with no other changes.
      These definitions were interchanged up to 10 Sept 1993.
      The wording was revised at this date.

   d. Re PP Field Code 2061 allocated below for lightning : 
      The SFERICS lightning count is expressed as an count per area per time.  
      The area is the area of a grid box e.g. on the verification grid. The time 
      is commonly one hour, but may be different depending on user chosen processing 
      options.

   e. Please Note : The range for Helicopter-Triggered lightning changed from 0-10 to 0-30 
      recently (03/07/12), thus, any data produced prior to 03/07/12 has a range of (0 - 10)

      The following documentation giving more detail about helicopter triggered lightning
      is available :-
      Investigation and prediction of helicopter triggered lightning over the North Sea.
      Wilkinson, J. M., H. Wells, P.R. Field, and P. Agnew Meteorol. Applications. 
      (early online release) http://onlinelibrary.wiley.com/doi/10.1002/met.1314/abstract

   f. The following online documentation for Lightning Index and Lightning Risk is available (internally
      at the Met Office):
      http://fcm9/projects/PostProc/wiki/PostProcDocConvection#Lightningindex01or10

   g. Updates to this pp/field codes document will propogate through to the `reference service <​http://reference.metoffice.gov.uk/um/fieldcode>`_.
 
