! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP.
!
! Variables for lbc addressing, sizes etc


MODULE lbc_mod
USE field_types, ONLY: &
    nfld_max
USE UM_Parparams, ONLY: &
    Nhalo_max
USE rimtypes, ONLY: &
    nrima_max
USE missing_data_mod, ONLY: &
    imdi

IMPLICIT NONE

PRIVATE :: imdi

! data structure sizes for atmos & ocean boundary file control
! routines
INTEGER :: rimwidtha(nrima_max)
INTEGER :: rimwidtho

! size of atmos lbc for given field type, halo type and rimwidth
! type
INTEGER :: lenrima(nfld_max,nhalo_max,nrima_max)

! size of atmos lbc on a given processor for given field type, halo type 
! and rimwidth type
! allocate at init, dims:(nfld_max,nhalo_max,nrima_max,0:nproc_max-1)
INTEGER, POINTER :: g_lenrima(:,:,:,:)=>NULL()

INTEGER :: nrim_timesa = imdi  ! Max no of timelevels in rim flds

! size of given side (pnorth,peast,psouth and pwest), field type,
! halo type and rimwidth type
INTEGER :: lbc_sizea(4,nfld_max,nhalo_max,nrima_max)

! start of a given side within the lbc
INTEGER :: lbc_starta(4,nfld_max,nhalo_max,nrima_max)

! start of a given side within the lbc on a given processor
! allocate at init, dims:(4,nfld_max,nhalo_max,nrima_max,0:maxproc-1)
INTEGER, POINTER :: g_lbc_starta(:,:,:,:,:)=>NULL()

! global (within the file) information

! size of atmos lbc on disk for given field type, halo type and
! rimwidth type
INTEGER :: global_lenrima(nfld_max,nhalo_max,nrima_max)

! size of given side, field type and halo type
INTEGER :: global_lbc_sizea(4,nfld_max,nhalo_max,nrima_max)

! start of a given side within the lbc
INTEGER :: global_lbc_starta(4,nfld_max,nhalo_max,nrima_max)

! variables that may be needed for vn5.2 but have not yet been
! dealt with at vn5.1
INTEGER :: rimfldsa
INTEGER :: boundflds
INTEGER :: global_lenrimdata_a
INTEGER :: lenrimdata_a
INTEGER :: rim_lookupsa
INTEGER :: bound_lookupsa

END MODULE lbc_mod
