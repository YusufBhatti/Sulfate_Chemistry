#if !defined(C_PIO_TIMER_H)
#define C_PIO_TIMER_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* Start C_IO_TIMER                                                 */
/*                                                                  */
/* Description:                                                     */
/*                                                                  */
/* Global constants I/O timer routines                              */
/*                                                                  */
/* Information:                                                     */
/*                                                                  */
/* Defines the variables required for IO timing + stats             */
/* Done globally as a number of functions will need them            */
/*                                                                  */
/*                                                                  */

/* DEPENDS ON: pio_io_timer.o */

#include <sys/time.h>
#include <unistd.h>
#include "c_fort2c_prototypes.h"

/* name mangling */

#if defined(C_LOW)
#define io_output_inter_timings io_output_inter_timings
#define io_timing_init io_timing_init
#elif defined(C_LOW_U)
#define io_output_inter_timings io_output_inter_timings_
#define io_timing_init io_timing_init_
#else
#define io_output_inter_timings IO_OUTPUT_INTER_TIMINGS
#define io_timing_init IO_TIMING_INIT
#endif

/* Constants */

#define IO_ITER 10

#define IO_B_IN  0 /* Buffin */
#define IO_B_OUT 1 /* Buffout */
#define IO_F_OP  2 /* File Open */
#define IO_F_CL  3 /* File Close */
#define IO_F_SK  4 /* File Seek */

#define IO_NUM_TIMERS 5

/* Prototypes */

extern void io_update_timer         (struct timeval,
                                     struct timeval,
                                     u_integer,
                                     integer);

extern void io_total_timings        (void);

extern void io_output_inter_timings (void);

extern void io_timing_init          (void);

/* Variables */

extern u_integer io_timer_active;

/* End C_IO_TIMER                                                   */

#endif
