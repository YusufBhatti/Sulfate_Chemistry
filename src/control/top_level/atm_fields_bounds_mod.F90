! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  defines a type to address array bounds for UM arrays
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

MODULE atm_fields_bounds_mod

USE model_domain_mod, ONLY: model_type, mt_single_column

IMPLICIT NONE
!
! Description:
!  switchable variables for loop bounds and array size declarations.
!  now only EG choices included.
!  v-at-the-poles and u-fields starting at 0
!  Default  initiation set to a value which will hopefully
!  crash any attempt of using these module variables uninitialised
!
! Method:
!
!
! Code description:
!   Language: Fortran 90.
!   This code is written to programming standard UMDP3 vn 8.1.


TYPE array_dims
  INTEGER :: i_start =-HUGE(INT(1))
  INTEGER :: i_end   = HUGE(INT(1))
  INTEGER :: i_len   = HUGE(INT(1))
  INTEGER :: j_start =-HUGE(INT(1))
  INTEGER :: j_end   = HUGE(INT(1))
  INTEGER :: j_len   = HUGE(INT(1))
  INTEGER :: k_start =-HUGE(INT(1))
  INTEGER :: k_end   = HUGE(INT(1))
  INTEGER :: k_len   = HUGE(INT(1))
  INTEGER :: halo_i  = HUGE(INT(1))
  INTEGER :: halo_j  = HUGE(INT(1))
END TYPE array_dims


! arrays for u,v,w,t,p fields without halo, small halo and large
! halo

TYPE (array_dims), SAVE                                                &
                  :: udims, vdims, wdims, tdims, pdims,                &
                     udims_s, vdims_s, wdims_s, tdims_s, pdims_s,      &
                     udims_l, vdims_l, wdims_l, tdims_l, pdims_l


! further arrays to cope with different number of vertical levels
TYPE (array_dims), SAVE                                                &
                        :: rdims2, o3dims2,                            &
                           o3dims, oneddims

! arrays used in stochastic physics (atm section=35) to define a
! smoothing mask with halo size set in stph_setup
TYPE (array_dims), SAVE :: stphdims_l

! Dimensions for SCM specific arrays
INTEGER :: ScmRowLen
INTEGER :: ScmRow

CONTAINS


SUBROUTINE atm_fields_bounds_init(offx,offy,halo_i,halo_j,             &
                  row_length,rows,n_rows,                              &
                  tr_levels_opt,bl_levels_opt,ozone_levels_opt)

USE nlsizes_namelist_mod, ONLY: model_levels

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                 &
       offx,offy,halo_i,halo_j,row_length,rows,n_rows

INTEGER, INTENT(IN), OPTIONAL ::                                       &
       bl_levels_opt,tr_levels_opt,ozone_levels_opt
INTEGER ::                                                             &
       bl_levels,tr_levels,ozone_levels


! Initialise to model levels
bl_levels    = model_levels
tr_levels    = model_levels
ozone_levels = model_levels

! If options are present then set the required levels.
IF (PRESENT(bl_levels_opt    )) bl_levels    = bl_levels_opt
IF (PRESENT(tr_levels_opt    )) tr_levels    = tr_levels_opt
IF (PRESENT(ozone_levels_opt )) ozone_levels = ozone_levels_opt

! Initialise SCM dimensions
! By default this is 1, though could be set in rose stem if
! scm evolved to have multiple columns
ScmRowLen = 1
ScmRow    = 1

! dimensions of u-field:
udims%i_start           = 0
udims%i_end             = row_length-1
  
IF (model_type == mt_single_column) THEN
  udims%i_start         = 1
  udims%i_end           = row_length
END IF 

udims%i_len             = udims%i_end - udims%i_start + 1
udims%j_start           = 1
udims%j_end             = rows
udims%j_len             = udims%j_end - udims%j_start + 1
udims%k_start           = 1
udims%k_end             = model_levels
udims%k_len             = udims%k_end - udims%k_start + 1
udims%halo_i            = 0
udims%halo_j            = 0

! dimensions of v-field:
vdims%i_start           = 1
vdims%i_end             = row_length
vdims%i_len             = vdims%i_end - vdims%i_start + 1
vdims%j_start           = 0
vdims%j_end             = n_rows-1

IF (model_type == mt_single_column) THEN
  vdims%j_start         = 1
  vdims%j_end           = n_rows
END IF 

vdims%j_len             = vdims%j_end - vdims%j_start + 1
vdims%k_start           = 1
vdims%k_end             = model_levels
vdims%k_len             = vdims%k_end - vdims%k_start + 1
vdims%halo_i            = 0
vdims%halo_j            = 0

! dimensions of w-field:
wdims%i_start           = 1
wdims%i_end             = row_length
wdims%i_len             = wdims%i_end - wdims%i_start + 1
wdims%j_start           = 1
wdims%j_end             = rows
wdims%j_len             = wdims%j_end - wdims%j_start + 1
wdims%k_start           = 0
wdims%k_end             = model_levels
wdims%k_len             = wdims%k_end - wdims%k_start + 1
wdims%halo_i            = 0
wdims%halo_j            = 0

! dimensions of theta-field:
tdims                   = wdims

! dimensions of pressure-field:
pdims%i_start           = 1
pdims%i_end             = row_length
pdims%i_len             = pdims%i_end - pdims%i_start + 1
pdims%j_start           = 1
pdims%j_end             = rows
pdims%j_len             = pdims%j_end - pdims%j_start + 1
pdims%k_start           = 1
pdims%k_end             = model_levels
pdims%k_len             = pdims%k_end - pdims%k_start + 1
pdims%halo_i            = 0
pdims%halo_j            = 0

! small halos

! dimensions of u-field (small halos):
udims_s%i_start         = -offx
udims_s%i_end           = row_length-1+offx

IF (model_type == mt_single_column) THEN
  udims_s%i_start       = 1-offx
  udims_s%i_end         = row_length+offx
END IF 

udims_s%i_len           = udims_s%i_end - udims_s%i_start + 1
udims_s%j_start         = 1-offy
udims_s%j_end           = rows+offy
udims_s%j_len           = udims_s%j_end - udims_s%j_start + 1
udims_s%k_start         = 1
udims_s%k_end           = model_levels
udims_s%k_len           = udims_s%k_end - udims_s%k_start + 1
udims_s%halo_i            = offx
udims_s%halo_j            = offy

! dimensions of v-field (small halos):
vdims_s%i_start         = 1-offx
vdims_s%i_end           = row_length+offx
vdims_s%i_len           = vdims_s%i_end - vdims_s%i_start + 1
vdims_s%j_start         = -offy
vdims_s%j_end           = n_rows-1+offy

IF (model_type == mt_single_column) THEN
  vdims_s%j_start       = 1-offy
  vdims_s%j_end         = n_rows+offy
END IF 

vdims_s%j_len           = vdims_s%j_end - vdims_s%j_start + 1
vdims_s%k_start         = vdims%k_start
vdims_s%k_end           = vdims%k_end
vdims_s%k_len           = vdims_s%k_end - vdims_s%k_start + 1
vdims_s%halo_i            = offx
vdims_s%halo_j            = offy

! dimensions of w-field (small halos):
wdims_s%i_start         = 1-offx
wdims_s%i_end           = row_length+offx
wdims_s%i_len           = wdims_s%i_end - wdims_s%i_start + 1
wdims_s%j_start         = 1-offy
wdims_s%j_end           = rows+offy
wdims_s%j_len           = wdims_s%j_end - wdims_s%j_start + 1
wdims_s%k_start         = 0
wdims_s%k_end           = model_levels
wdims_s%k_len           = wdims_s%k_end - wdims_s%k_start + 1
wdims_s%halo_i          = offx
wdims_s%halo_j          = offy

! dimensions of theta-field (small halos):
tdims_s                 = wdims_s

! dimensions of pressure-field (small halos):
pdims_s%i_start         = 1-offx
pdims_s%i_end           = row_length+offx
pdims_s%i_len           = pdims_s%i_end - pdims_s%i_start + 1
pdims_s%j_start         = 1-offy
pdims_s%j_end           = rows+offy
pdims_s%j_len           = pdims_s%j_end - pdims_s%j_start + 1
pdims_s%k_start         = 1
pdims_s%k_end           = model_levels
pdims_s%k_len           = pdims_s%k_end - pdims_s%k_start + 1
pdims_s%halo_i          = offx
pdims_s%halo_j          = offy

! large halos
! dimensions of u-field (large halos):
udims_l%i_start         = -halo_i
udims_l%i_end           = row_length-1+halo_i

IF (model_type == mt_single_column) THEN
  udims_l%i_start       = 1-halo_i
  udims_l%i_end         = row_length+halo_i
END IF 

udims_l%i_len           = udims_l%i_end - udims_l%i_start + 1
udims_l%j_start         = 1-halo_j
udims_l%j_end           = rows+halo_j
udims_l%j_len           = udims_l%j_end - udims_l%j_start + 1
udims_l%k_start         = 1
udims_l%k_end           = model_levels
udims_l%k_len           = udims_l%k_end - udims_l%k_start + 1
udims_l%halo_i          = halo_i
udims_l%halo_j          = halo_j

! dimensions of v-field (large halos):
vdims_l%i_start         = 1-halo_i
vdims_l%i_end           = row_length+halo_i
vdims_l%i_len           = vdims_l%i_end - vdims_l%i_start + 1
vdims_l%j_start         = -halo_j
vdims_l%j_end           = n_rows-1+halo_j

IF (model_type == mt_single_column) THEN
  vdims_l%j_start       = 1-halo_j
  vdims_l%j_end         = n_rows+halo_j
END IF 

vdims_l%j_len           = vdims_l%j_end - vdims_l%j_start + 1
vdims_l%k_start         = vdims%k_start
vdims_l%k_end           = vdims%k_end
vdims_l%k_len           = vdims_l%k_end - vdims_l%k_start + 1
vdims_l%halo_i          = halo_i
vdims_l%halo_j          = halo_j

! dimensions of w-field (large halos):
wdims_l%i_start         = 1-halo_i
wdims_l%i_end           = row_length+halo_i
wdims_l%i_len           = wdims_l%i_end - wdims_l%i_start + 1
wdims_l%j_start         = 1-halo_j
wdims_l%j_end           = rows+halo_j
wdims_l%j_len           = wdims_l%j_end - wdims_l%j_start + 1
wdims_l%k_start         = wdims%k_start
wdims_l%k_end           = wdims%k_end
wdims_l%k_len           = wdims_l%k_end - wdims_l%k_start + 1
wdims_l%halo_i          = halo_i
wdims_l%halo_j          = halo_j

! dimensions of theta-field (large halos):
tdims_l                 = wdims_l

! dimensions of pressure-field (large halos):
pdims_l%i_start         = 1-halo_i
pdims_l%i_end           = row_length+halo_i
pdims_l%i_len           = pdims_l%i_end - pdims_l%i_start + 1
pdims_l%j_start         = 1-halo_j
pdims_l%j_end           = rows+halo_j
pdims_l%j_len           = pdims_l%j_end - pdims_l%j_start + 1
pdims_l%k_start         = 1
pdims_l%k_end           = model_levels
pdims_l%k_len           = pdims_l%k_end - pdims_l%k_start + 1
pdims_l%halo_i          = halo_i
pdims_l%halo_j          = halo_j

! other fields
! dimensions of rho-field :
rdims2                  = pdims
rdims2%k_start          = 0
rdims2%k_end            = model_levels+1
rdims2%k_len            = rdims2%k_end - rdims2%k_start + 1

! dimensions of O3-array
o3dims2                 = tdims

! Ozone dimensions (1d):
o3dims%i_start          = 1
o3dims%i_end            = row_length*rows*(model_levels+1)
o3dims%i_len            = o3dims%i_end - o3dims%i_start + 1
o3dims%j_start          = 0
o3dims%j_end            = 0
o3dims%j_len            = o3dims%j_end - o3dims%j_start + 1
o3dims%k_start          = 0
o3dims%k_end            = 0
o3dims%k_len            = o3dims%k_end - o3dims%k_start + 1

oneddims%i_start          = 1
oneddims%i_end            = row_length*rows*model_levels
oneddims%i_len            = oneddims%i_end - oneddims%i_start + 1
oneddims%j_start          = 0
oneddims%j_end            = 0
oneddims%j_len            = oneddims%j_end - oneddims%j_start + 1
oneddims%k_start          = 0
oneddims%k_end            = 0
oneddims%k_len            = oneddims%k_end - oneddims%k_start + 1

END SUBROUTINE atm_fields_bounds_init
END MODULE atm_fields_bounds_mod

