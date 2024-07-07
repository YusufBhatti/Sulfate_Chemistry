! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Decription:
! given an integer's data_start,data_end, splits this range into equal 
! chunk_sizes for each thread. If the range is not divisible by the number
! of threads, chunk_size is increased by one for all the threads, starting 
! with the last and counting backwards. Also returns start and end as indices
! into a shared array. Procedures can be made thread safe by calling this 
! routine to determine the worksharing pattern. They can then be called from 
! a parallel region with private variables.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
! 
MODULE compute_chunk_size_mod

  IMPLICIT NONE
CONTAINS
  
  SUBROUTINE compute_chunk_size(data_start,data_end,iter_start,iter_end)
!$    USE OMP_lib 
! only use omp_lib if compiling with openMP 
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: data_start!start of range to be split across threads
    INTEGER, INTENT(IN)  :: data_end  !end of range to be split across threads
    INTEGER, INTENT(OUT) :: iter_start ! start for this thread
    INTEGER, INTENT(OUT) :: iter_end ! end for this thread
    ! local variables
    INTEGER :: thread_id, num_threads, extra
    INTEGER :: chunk_size
    INTEGER :: data_range  ! range of values from start to end

    ! outside of parallel region num_threads will be 1 if compiled without omp
    ! num_threads will be 1
    num_threads = 1
!$     num_threads=omp_get_num_threads()

    thread_id = 0
!$    thread_id=omp_get_thread_num()

    data_range = (data_end - data_start) + 1
    chunk_size = data_range/num_threads
    extra = data_range - chunk_size*num_threads
       
    iter_start = chunk_size*thread_id + data_start
    iter_end = iter_start + chunk_size - 1

    ! account for data_range not divisible by num_threads
    IF( (thread_id >= (num_threads-extra)) .AND. (extra /= 0) ) THEN
       chunk_size = chunk_size + 1
       iter_start = iter_start + thread_id - (num_threads-extra)
       iter_end = iter_end + thread_id - (num_threads-extra) + 1
    END IF
    
    RETURN
  END SUBROUTINE compute_chunk_size
END MODULE compute_chunk_size_mod
