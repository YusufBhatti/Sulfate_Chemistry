! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! This module provide the tags used by the various point-to-point
! MPI/MPL messages in the UM

MODULE tags_params

IMPLICIT NONE

! -----------------------------------------------------------------!
! REMINDER:                                                        !
! The same tag can be used safely with different MPI communicators !
! -----------------------------------------------------------------!

! Info         : Tags used for halo exchange
! Communicator : halos_comm
! Range        : 1:20
! Note         : 

INTEGER, PARAMETER :: tag_to_e     = 1       
INTEGER, PARAMETER :: tag_to_w     = 2       
INTEGER, PARAMETER :: tag_from_e   = tag_to_w
INTEGER, PARAMETER :: tag_from_w   = tag_to_e
INTEGER, PARAMETER :: tag_to_s     = 3           ! not polar  
INTEGER, PARAMETER :: tag_to_n     = 4           ! not polar
INTEGER, PARAMETER :: tag_from_s   = tag_to_n    ! not polar
INTEGER, PARAMETER :: tag_from_n   = tag_to_s    ! not polar
INTEGER, PARAMETER :: tag_to_s_p   = 13          ! polar  
INTEGER, PARAMETER :: tag_to_n_p   = 14          ! polar
INTEGER, PARAMETER :: tag_from_s_p = tag_to_s_p  ! polar
INTEGER, PARAMETER :: tag_from_n_p = tag_to_n_p  ! polar
INTEGER, PARAMETER :: tag_to_ne    = 5         
INTEGER, PARAMETER :: tag_to_se    = 6         
INTEGER, PARAMETER :: tag_to_sw    = 7         
INTEGER, PARAMETER :: tag_to_nw    = 8         
INTEGER, PARAMETER :: tag_from_ne  = tag_to_sw   
INTEGER, PARAMETER :: tag_from_se  = tag_to_nw   
INTEGER, PARAMETER :: tag_from_sw  = tag_to_ne   
INTEGER, PARAMETER :: tag_from_nw  = tag_to_se   
INTEGER, PARAMETER :: tag_to_ne_p  = 15          ! polar
INTEGER, PARAMETER :: tag_to_se_p  = 16          ! polar
INTEGER, PARAMETER :: tag_to_sw_p  = 17          ! polar
INTEGER, PARAMETER :: tag_to_nw_p  = 18          ! polar
INTEGER, PARAMETER :: tag_from_ne_p= tag_to_nw_p ! polar
INTEGER, PARAMETER :: tag_from_nw_p= tag_to_ne_p ! polar  
INTEGER, PARAMETER :: tag_from_se_p= tag_to_sw_p ! polar
INTEGER, PARAMETER :: tag_from_sw_p= tag_to_se_p ! polar  


! Info         : Tags used for non-blocking halo exchange
! Communicator : halos_comm
! Range        : 101:32020
! Note         :  
!
! The tags for the non blocking halo exchange are composed by:
! The halo tag (1:20) + non blocking tag base (100) * non blocking tag index.
! Each non-blocking swap bounds is identified by an unique tag index.
!

! The base value to be added to the tags used by the swap_bounds.
INTEGER, PARAMETER :: nb_tag_base = 100

! The maximum value allowed to the global tag_index.
! The tag index is multiplied to the tag base
! MPI standard requires that the maximum tag value is at least 32767
INTEGER, PARAMETER :: nb_max_tag_index  = 320


END MODULE tags_params
