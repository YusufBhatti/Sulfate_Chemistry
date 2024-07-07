! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Identify required climatology fields and fields to pass to RADAER
!           within the GLOMAP_CLIM Section
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE glomap_clim_identify_fields_mod

IMPLICIT NONE

! -----------------------------------------------------------------------------*
!                                                                              *
! Procedure:                                                                   *
!   1) Using settings of mode and component (ukca_mode_setup), climatology     *
!       aerosol fields are identified as being required by RADAER.             *
!                                                                              *
!   2) Using settings of mode and component (ukca_mode_setup), other input     *
!       fields are identified as being required by RADAER.                     *
!                                                                              *
! -----------------------------------------------------------------------------*

PUBLIC
SAVE

TYPE codetype
  INTEGER :: section        ! section code
  INTEGER :: item           ! item code
  LOGICAL :: put_stash      ! item has to be written to stash
  LOGICAL :: required       ! t/f
ENDTYPE codetype

! List of fields - either climatology aerosol fields or fields to pass to RADAER
TYPE(codetype), ALLOCATABLE :: aerofields(:)

! Length of aerofields array
INTEGER, PARAMETER :: n_stored_items = 150 

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName =                           &
                                              'GLOMAP_CLIM_IDENTIFY_FIELDS_MOD'

CONTAINS

! #############################################################################

SUBROUTINE glomap_clim_identify_fields

USE missing_data_mod,        ONLY: &
    imdi

USE parkind1,                ONLY: &
    jprb,                          &
    jpim

USE ukca_mode_setup,         ONLY: &
    mode,                          &
    component,                     &
    cp_su,                         &
    cp_bc,                         &
    cp_oc,                         &
    cp_cl,                         &
    cp_du,                         &
    cp_so,                         &
    mode_nuc_sol,                  &
    mode_ait_sol,                  &
    mode_acc_sol,                  &
    mode_cor_sol,                  &
    mode_ait_insol,                &
    mode_acc_insol,                &
    mode_cor_insol

USE um_stashcode_mod,        ONLY: &
    stashcode_glomap_clim_sec,     &
    stashcode_gc_nd_nuc_sol,       &
    stashcode_gc_nuc_sol_su,       &
    stashcode_gc_nuc_sol_oc,       &
    stashcode_gc_nd_ait_sol,       &
    stashcode_gc_ait_sol_su,       &
    stashcode_gc_ait_sol_bc,       &
    stashcode_gc_ait_sol_oc,       &
    stashcode_gc_nd_acc_sol,       &
    stashcode_gc_acc_sol_su,       &
    stashcode_gc_acc_sol_bc,       &
    stashcode_gc_acc_sol_oc,       &
    stashcode_gc_acc_sol_ss,       &
    stashcode_gc_acc_sol_du,       &
    stashcode_gc_nd_cor_sol,       &
    stashcode_gc_cor_sol_su,       &
    stashcode_gc_cor_sol_bc,       &
    stashcode_gc_cor_sol_oc,       &
    stashcode_gc_cor_sol_ss,       &
    stashcode_gc_cor_sol_du,       &
    stashcode_gc_nd_ait_ins,       &
    stashcode_gc_ait_ins_bc,       &
    stashcode_gc_ait_ins_oc,       &
    stashcode_gc_nd_acc_ins,       &
    stashcode_gc_acc_ins_du,       &
    stashcode_gc_nd_cor_ins,       &
    stashcode_gc_cor_ins_du,       &
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

INTEGER, PARAMETER :: msect = 1000*stashcode_glomap_clim_sec  ! 1000 * section
INTEGER :: n_io_fields_p         ! No of aerosol fields to be read in or output
INTEGER :: n_io_fields_d         ! No of aerosol fields to be read in or output
INTEGER :: j                     ! Counter

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'GLOMAP_CLIM_IDENTIFY_FIELDS'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! List of Section 54 items to be read from D1 array / copied to stash
IF (.NOT. ALLOCATED(aerofields)) ALLOCATE(aerofields(n_stored_items))

! Do not include halo sizes here, only the leading dimensions are required
j = 1
n_io_fields_p = 26

! Set logical %required to false for GLOMAP-mode climatology aerosol fields,
!  then set below depending on the selected modes and components

aerofields(j:j+n_io_fields_p)%section    = stashcode_glomap_clim_sec
aerofields(:)%item                       = imdi
aerofields(:)%put_stash                  = .FALSE.
aerofields(:)%required                   = .FALSE.

! Aerosol number and mass mixing ratios, turn on as per modes and components

aerofields(j+0)%item=stashcode_gc_nd_nuc_sol - msect  ! Nucleation sol number
IF (mode(mode_nuc_sol))            THEN
  aerofields(j+0)%required  =  .TRUE.
END IF

aerofields(j+1)%item=stashcode_gc_nuc_sol_su - msect  ! Nucleation sol mode SU
IF (component(mode_nuc_sol,cp_su)) THEN
  aerofields(j+1)%required  =  .TRUE.
END IF

aerofields(j+25)%item=stashcode_gc_nuc_sol_oc - msect ! Nucleation sol mode OC
IF (component(mode_nuc_sol,cp_oc)) THEN
  aerofields(j+25)%required =  .TRUE.
END IF

aerofields(j+2)%item=stashcode_gc_nd_ait_sol - msect  ! Aitken sol number
IF (mode(mode_ait_sol))            THEN
  aerofields(j+2)%required  =  .TRUE.
END IF

aerofields(j+3)%item=stashcode_gc_ait_sol_su - msect  ! Aitken sol mode SU
IF (component(mode_ait_sol,cp_su)) THEN
  aerofields(j+3)%required  =  .TRUE.
END IF

aerofields(j+4)%item=stashcode_gc_ait_sol_bc - msect  ! Aitken sol mode BC
IF (component(mode_ait_sol,cp_bc)) THEN
  aerofields(j+4)%required  =  .TRUE.
END IF

aerofields(j+5)%item=stashcode_gc_ait_sol_oc - msect  ! Aitken sol mode OC
IF (component(mode_ait_sol,cp_oc)) THEN
  aerofields(j+5)%required  =  .TRUE.
END IF

aerofields(j+6)%item=stashcode_gc_nd_acc_sol - msect  ! Accumulation sol number
IF (mode(mode_acc_sol))            THEN
  aerofields(j+6)%required  =  .TRUE.
END IF

aerofields(j+7)%item=stashcode_gc_acc_sol_su - msect  ! Accumulation sol mode SU
IF (component(mode_acc_sol,cp_su)) THEN
  aerofields(j+7)%required  =  .TRUE.
END IF

aerofields(j+8)%item=stashcode_gc_acc_sol_bc - msect  ! Accumulation sol mode BC
IF (component(mode_acc_sol,cp_bc)) THEN
  aerofields(j+8)%required  =  .TRUE.
END IF

aerofields(j+9)%item=stashcode_gc_acc_sol_oc - msect  ! Accumulation sol mode OC
IF (component(mode_acc_sol,cp_oc)) THEN
  aerofields(j+9)%required =   .TRUE.
END IF

aerofields(j+10)%item=stashcode_gc_acc_sol_ss - msect ! Accumulation sol mode SS
IF (component(mode_acc_sol,cp_cl)) THEN
  aerofields(j+10)%required =  .TRUE.
END IF

aerofields(j+11)%item=stashcode_gc_acc_sol_du - msect ! Accumulation sol mode DU
IF (component(mode_acc_sol,cp_du)) THEN
  aerofields(j+11)%required =  .TRUE.
END IF

aerofields(j+12)%item=stashcode_gc_nd_cor_sol - msect ! Coarse sol number
IF (mode(mode_cor_sol))            THEN
  aerofields(j+12)%required =  .TRUE.
END IF

aerofields(j+13)%item=stashcode_gc_cor_sol_su - msect ! Coarse sol mode SU
IF (component(mode_cor_sol,cp_su)) THEN
  aerofields(j+13)%required =  .TRUE.
END IF

aerofields(j+14)%item=stashcode_gc_cor_sol_bc - msect ! Coarse sol mode BC
IF (component(mode_cor_sol,cp_bc)) THEN
  aerofields(j+14)%required =  .TRUE.
END IF

aerofields(j+15)%item=stashcode_gc_cor_sol_oc - msect ! Coarse sol mode OC
IF (component(mode_cor_sol,cp_oc)) THEN
  aerofields(j+15)%required =  .TRUE.
END IF

aerofields(j+16)%item=stashcode_gc_cor_sol_ss - msect ! Coarse sol mode SS
IF (component(mode_cor_sol,cp_cl)) THEN
  aerofields(j+16)%required =  .TRUE.
END IF

aerofields(j+17)%item=stashcode_gc_cor_sol_du - msect ! Coarse sol mode DU
IF (component(mode_cor_sol,cp_du)) THEN
  aerofields(j+17)%required =  .TRUE.
END IF

aerofields(j+18)%item=stashcode_gc_nd_ait_ins - msect ! Aitken ins number
IF (mode(mode_ait_insol))            THEN
  aerofields(j+18)%required =  .TRUE.
END IF

aerofields(j+19)%item=stashcode_gc_ait_ins_bc - msect ! Aitken ins mode BC
IF (component(mode_ait_insol,cp_bc)) THEN
  aerofields(j+19)%required =  .TRUE.
END IF

aerofields(j+20)%item=stashcode_gc_ait_ins_oc - msect ! Aitken ins mode OC
IF (component(mode_ait_insol,cp_oc)) THEN
  aerofields(j+20)%required =  .TRUE.
END IF

aerofields(j+21)%item=stashcode_gc_nd_acc_ins - msect ! Accumulation ins number
IF (mode(mode_acc_insol))            THEN
  aerofields(j+21)%required =  .TRUE.
END IF

aerofields(j+22)%item=stashcode_gc_acc_ins_du - msect ! Accumulation ins mode DU
IF (component(mode_acc_insol,cp_du)) THEN
  aerofields(j+22)%required =  .TRUE.
END IF

aerofields(j+23)%item=stashcode_gc_nd_cor_ins - msect ! Coarse ins number
IF (mode(mode_cor_insol))            THEN
  aerofields(j+23)%required =  .TRUE.
END IF

aerofields(j+24)%item=stashcode_gc_cor_ins_du - msect ! Coarse ins DU
IF (component(mode_cor_insol,cp_du)) THEN
  aerofields(j+24)%required =  .TRUE.
END IF
     
! Section 54 fields passed to RADAER
!  Item numbers correspond to those in routine: ukca_radaer_init
!  Set %put_stash where required as indicated by logical arrays mode, component

j = 99
n_io_fields_d = 38
aerofields(j+1:j+n_io_fields_d)%required   = .FALSE.
aerofields(j+1:j+n_io_fields_d)%put_stash  = .FALSE.
aerofields(j+1:j+n_io_fields_d)%section    = stashcode_glomap_clim_sec

! dry diameter
aerofields(j+1)%item = stashcode_gc_dryd_ait_sol - msect ! ait_sol
IF (mode(mode_ait_sol))    aerofields(j+1)%put_stash  =  .TRUE.

aerofields(j+2)%item = stashcode_gc_dryd_acc_sol - msect ! acc_sol
IF (mode(mode_acc_sol))    aerofields(j+2)%put_stash  =  .TRUE.

aerofields(j+3)%item = stashcode_gc_dryd_cor_sol - msect ! cor_sol
IF (mode(mode_cor_sol))    aerofields(j+3)%put_stash  =  .TRUE.

aerofields(j+4)%item = stashcode_gc_dryd_ait_ins - msect ! ait_ins
IF (mode(mode_ait_insol))  aerofields(j+4)%put_stash  =  .TRUE.

aerofields(j+5)%item = stashcode_gc_dryd_acc_ins - msect ! acc_ins
IF (mode(mode_acc_insol))  aerofields(j+5)%put_stash  =  .TRUE.

aerofields(j+6)%item = stashcode_gc_dryd_cor_ins - msect ! cor_ins
IF (mode(mode_cor_insol))  aerofields(j+6)%put_stash  =  .TRUE.

! wet radius
aerofields(j+7)%item = stashcode_gc_wetd_ait_sol - msect ! ait_sol
IF (mode(mode_ait_sol))    aerofields(j+7)%put_stash  =  .TRUE.

aerofields(j+8)%item = stashcode_gc_wetd_acc_sol - msect ! acc_sol
IF (mode(mode_acc_sol))    aerofields(j+8)%put_stash  =  .TRUE.

aerofields(j+9)%item = stashcode_gc_wetd_cor_sol - msect ! cor_sol
IF (mode(mode_cor_sol))    aerofields(j+9)%put_stash  =  .TRUE.

! aerosol density
aerofields(j+10)%item = stashcode_gc_rho_ait_sol - msect ! ait_sol
IF (mode(mode_ait_sol))    aerofields(j+10)%put_stash =  .TRUE.

aerofields(j+11)%item = stashcode_gc_rho_acc_sol - msect ! acc_sol
IF (mode(mode_acc_sol))    aerofields(j+11)%put_stash =  .TRUE.

aerofields(j+12)%item = stashcode_gc_rho_cor_sol - msect ! cor_sol
IF (mode(mode_cor_sol))    aerofields(j+12)%put_stash =  .TRUE.

aerofields(j+13)%item = stashcode_gc_rho_ait_ins - msect ! ait_ins
IF (mode(mode_ait_insol))  aerofields(j+13)%put_stash =  .TRUE.

aerofields(j+14)%item = stashcode_gc_rho_acc_ins - msect ! acc_ins
IF (mode(mode_acc_insol))  aerofields(j+14)%put_stash =  .TRUE.

aerofields(j+15)%item = stashcode_gc_rho_cor_ins - msect ! cor_ins
IF (mode(mode_cor_insol))  aerofields(j+15)%put_stash =  .TRUE.

! partial volume aitken soluble
aerofields(j+16)%item = stashcode_gc_pvol_ait_su_sol  - msect ! SU
IF (mode(mode_ait_sol) .AND. component(mode_ait_sol,cp_su))     THEN
  aerofields(j+16)%put_stash = .TRUE.
END IF

aerofields(j+17)%item = stashcode_gc_pvol_ait_bc_sol  - msect ! BC
IF (mode(mode_ait_sol) .AND. component(mode_ait_sol,cp_bc))     THEN
  aerofields(j+17)%put_stash = .TRUE.
END IF

aerofields(j+18)%item = stashcode_gc_pvol_ait_oc_sol  - msect ! OC
IF (mode(mode_ait_sol) .AND. component(mode_ait_sol,cp_oc))     THEN
  aerofields(j+18)%put_stash = .TRUE.
END IF

aerofields(j+19)%item = stashcode_gc_pvol_ait_so_sol  - msect ! SO
IF (mode(mode_ait_sol) .AND. component(mode_ait_sol,cp_so))     THEN
  aerofields(j+18)%put_stash = .TRUE.
END IF

aerofields(j+20)%item = stashcode_gc_pvol_ait_h2o_sol - msect ! H2O
IF (mode(mode_ait_sol))   aerofields(j+20)%put_stash =  .TRUE.

! partial volume accumulation soluble
aerofields(j+21)%item = stashcode_gc_pvol_acc_su_sol  - msect ! SU
IF (mode(mode_acc_sol) .AND. component(mode_acc_sol,cp_su))     THEN
  aerofields(j+21)%put_stash = .TRUE.
END IF

aerofields(j+22)%item = stashcode_gc_pvol_acc_bc_sol  - msect ! BC
IF (mode(mode_acc_sol) .AND. component(mode_acc_sol,cp_bc))     THEN
  aerofields(j+22)%put_stash = .TRUE.
END IF

aerofields(j+23)%item = stashcode_gc_pvol_acc_oc_sol  - msect ! OC
IF (mode(mode_acc_sol) .AND. component(mode_acc_sol,cp_oc))     THEN
  aerofields(j+23)%put_stash = .TRUE.
END IF

aerofields(j+24)%item = stashcode_gc_pvol_acc_ss_sol  - msect ! SS
IF (mode(mode_acc_sol) .AND. component(mode_acc_sol,cp_cl))     THEN
  aerofields(j+24)%put_stash = .TRUE.
END IF

aerofields(j+25)%item = stashcode_gc_pvol_acc_du_sol  - msect ! DU
IF (mode(mode_acc_sol) .AND. component(mode_acc_sol,cp_du))     THEN
  aerofields(j+25)%put_stash = .TRUE.
END IF

aerofields(j+26)%item = stashcode_gc_pvol_acc_so_sol  - msect ! SO
IF (mode(mode_acc_sol) .AND. component(mode_acc_sol,cp_so))     THEN
  aerofields(j+26)%put_stash = .TRUE.
END IF

aerofields(j+27)%item = stashcode_gc_pvol_acc_h2o_sol - msect ! H2O
IF (mode(mode_acc_sol))   aerofields(j+27)%put_stash =  .TRUE.

! partial volume coarse soluble
aerofields(j+28)%item = stashcode_gc_pvol_cor_su_sol  - msect ! SU
IF (mode(mode_cor_sol) .AND. component(mode_cor_sol,cp_su))     THEN
  aerofields(j+28)%put_stash = .TRUE.
END IF

aerofields(j+29)%item = stashcode_gc_pvol_cor_bc_sol  - msect ! BC
IF (mode(mode_cor_sol) .AND. component(mode_cor_sol,cp_bc))     THEN
  aerofields(j+29)%put_stash = .TRUE.
END IF

aerofields(j+30)%item = stashcode_gc_pvol_cor_oc_sol  - msect ! OC
IF (mode(mode_cor_sol) .AND. component(mode_cor_sol,cp_oc))     THEN
  aerofields(j+30)%put_stash = .TRUE.
END IF

aerofields(j+31)%item = stashcode_gc_pvol_cor_ss_sol  - msect ! SS
IF (mode(mode_cor_sol) .AND. component(mode_cor_sol,cp_cl))     THEN
  aerofields(j+31)%put_stash = .TRUE.
END IF

aerofields(j+32)%item = stashcode_gc_pvol_cor_du_sol  - msect ! DU
IF (mode(mode_cor_sol) .AND. component(mode_cor_sol,cp_du))     THEN
  aerofields(j+32)%put_stash = .TRUE.
END IF

aerofields(j+33)%item = stashcode_gc_pvol_cor_so_sol  - msect ! SO
IF (mode(mode_cor_sol) .AND. component(mode_cor_sol,cp_so))     THEN
  aerofields(j+33)%put_stash = .TRUE.
END IF

aerofields(j+34)%item = stashcode_gc_pvol_cor_h2o_sol - msect ! H2O
IF (mode(mode_cor_sol))   aerofields(j+34)%put_stash =  .TRUE.

! partial volume aitken insoluble
aerofields(j+35)%item = stashcode_gc_pvol_ait_bc_ins  - msect ! BC
IF (mode(mode_ait_insol) .AND. component(mode_ait_insol,cp_bc)) THEN
  aerofields(j+35)%put_stash = .TRUE.
END IF

aerofields(j+36)%item = stashcode_gc_pvol_ait_oc_ins  - msect ! OC
IF (mode(mode_ait_insol) .AND. component(mode_ait_insol,cp_oc)) THEN
  aerofields(j+36)%put_stash = .TRUE.
END IF

! partial volume accumulation insoluble
aerofields(j+37)%item = stashcode_gc_pvol_acc_du_ins  - msect ! DU
IF (mode(mode_acc_insol) .AND. component(mode_acc_insol,cp_du)) THEN
  aerofields(j+37)%put_stash = .TRUE.
END IF

! partial volume coarse insoluble
aerofields(j+38)%item = stashcode_gc_pvol_cor_du_ins  - msect ! DU
IF (mode(mode_cor_insol) .AND. component(mode_cor_insol,cp_du)) THEN
  aerofields(j+38)%put_stash = .TRUE.
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE glomap_clim_identify_fields

END MODULE glomap_clim_identify_fields_mod
