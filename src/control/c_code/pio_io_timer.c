#if defined(C95_2A)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/**********************************************************************/
/* Description:                                                       */
/*                                                                    */
/* I/O Timing Routines                                                */
/*                                                                    */
/* Information:                                                       */
/*                                                                    */
/* Provides routines for I/O timing. The following routines are       */
/* available:                                                         */
/*                                                                    */
/*   io_eval_overhead     - Called from C.  This function attempts to */
/*                          establish how long a call to gettimeofday */
/*                          takes.                                    */
/*   io_update_timer      - Called from C. Wrapper to call the actual */
/*                          timer update function.                    */
/*   io_update_timer_data - Called from C. Function to update the     */
/*                          cumulative timers                         */
/*   io_report            - Called from C. Function to provide a      */
/*                          simple report for I/O timing.             */
/*   io_total_timings     - Called from Fortran.  Function to produce */
/*                          a simple report detailing I/O usage for   */
/*                          the complete program to date.             */
/*   io_output_inter_timings  - Called from Fortran. Function to      */
/*                              produce a simple report detailing the */
/*                              I/O usage since the last call of this */
/*                              function.                             */
/*   io_timing_init       - Called from Fortran.  Function to set the */
/*                          value of the io_timer_active flag (0 if   */
/*                          disabled; non-zero if active).            */
/*                                                                    */
/* For further details, see UMDP S5.                                  */
/*                                                                    */
/* Code Description:                                                  */
/*   Language: C                                                      */
/*   This code is written to UMDP3 v6 programming standards.          */
/*                                                                    */
/**********************************************************************/

/* Standard header files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Header files for Fortran to C interface */
#include "c_fort2c_prototypes.h"
#include "c_fort2c_types.h"
#include "pio_umprint.h"

/* Header files for IO Timing */
#include "c_pio_timer.h"

/* Global variables */
static u_integer io_calls[IO_NUM_TIMERS] = {0};  /* No. of calls to routine */
static u_integer io_bytes[IO_NUM_TIMERS] = {0};  /* Bytes transferred       */
static u_integer io_secs[IO_NUM_TIMERS]  = {0};  /* Seconds taken           */
static u_integer io_usecs[IO_NUM_TIMERS] = {0};  /* Microseconds taken      */

static u_integer io_calls_last[IO_NUM_TIMERS] = {0};  /* Previous values of */
static u_integer io_bytes_last[IO_NUM_TIMERS] = {0};  /* IO data for        */
static u_integer io_secs_last[IO_NUM_TIMERS]  = {0};  /* computing          */
static u_integer io_usecs_last[IO_NUM_TIMERS] = {0};  /* incrementals       */

static u_integer io_overhead       = 0;       /* timer overhead       */

u_integer io_timer_active   = 0;       /* main IO timer switch */


/* Prototypes for internally called functions */

static void io_eval_overhead     (void);

static void io_update_timer_data (struct timeval,
                                  struct timeval,
                                  u_integer,
                                  u_integer *,
                                  u_integer *,
                                  u_integer *,
                                  u_integer *);

static void io_report            (u_integer,
                                  u_integer,
                                  u_integer,
                                  u_integer,
                                  integer);

/* The following routines are for IO timing */
/****************************************************************/
/* The following function trys to make a guess at how long a    */
/* call to gettimeofday takes. It does this by making two calls */
/* immediately adjacent and taking the difference. This is      */
/* averaged over a number of calls to improve the figure given  */
/****************************************************************/
static void io_eval_overhead( void )
{
  int i;
  struct timeval first,
                 second;

  /* if overhead is non-zero it has already been calculated. */
  if ( io_overhead != 0 )
    return;

  /* Loop - calling gettimeofday the required no. of times */
  for (i=0; i<IO_ITER; i++)
  {
    gettimeofday( &first,  NULL );
    gettimeofday( &second, NULL );

    io_overhead += (u_integer)( (second.tv_sec - first.tv_sec) * 1000000
                         +((integer)second.tv_usec - (integer)first.tv_usec) );
  }

  io_overhead /= IO_ITER;
}


/****************************************************************/
/* The following function is a wrapper to call the actual timer */
/* update function for both total and last data.                */
/****************************************************************/
void io_update_timer( struct timeval first,
                      struct timeval second,
                      u_integer bytes_to_add,
                      integer op )
{
  io_update_timer_data( first, second, bytes_to_add,
                        &(io_calls[op]), &(io_bytes[op]),
                        &(io_secs[op]),  &(io_usecs[op]) );

  io_update_timer_data( first, second, bytes_to_add,
                        &(io_calls_last[op]), &(io_bytes_last[op]),
                        &(io_secs_last[op]),  &(io_usecs_last[op]) );
  return;
}

/****************************************************************/
/* The following function is used to update the cumulative      */
/* timers etc that are held.                                    */
/****************************************************************/
static void io_update_timer_data( struct timeval first,
                           struct timeval second,
                           u_integer bytes_to_add,
                           u_integer *calls,
                           u_integer *bytes,
                           u_integer *secs,
                           u_integer *usecs)
{
  /* Update the timers */
  (*secs)  += (u_integer)( second.tv_sec  - first.tv_sec  );

  /* Cater for roll over of second */
  if (first.tv_usec > second.tv_usec)
  {
    (*usecs) += 1000000;
    (*secs)--;
  }

  (*usecs) = (u_integer)((integer)(*usecs) + ((integer)second.tv_usec
                                                     - (integer)first.tv_usec));

  /* cater for the function call overhead */
  io_eval_overhead();
  (*usecs) -= io_overhead;

  if ( *usecs  > 1000000 )
  {
    (*usecs) -= 1000000;
    (*secs) ++;
  }

  /* increment bytes and calls */
  (*calls) ++;
  (*bytes) += bytes_to_add;
}

/****************************************************************/
/* The following function is a utility to write out a simple    */
/* report for the IO                                            */
/****************************************************************/
static void io_report(u_integer calls,
                      u_integer bytes,
                      u_integer secs,
                      u_integer usecs,
                      integer op )
{
  u_integer avg_secs;
  u_integer avg_usecs;
  u_integer avg_bytes;

  float time;
  float avg_time;
  float mbs;         /* Megabytes per second */
  float avg_mbs;     /* Average mbs */

  char message[MAX_OUTSTR];

  /* Calculate the averages we need for output */
  if ( calls == 0 )
  {
    avg_usecs = 0;
    avg_secs  = 0;
    avg_bytes = 0;
  }
  else
  {
    avg_usecs = ((1000000 * secs) + usecs) / calls;
    avg_secs = avg_usecs / 1000000;
    avg_usecs = avg_usecs % 1000000;

    avg_bytes = bytes/calls;
  }

  /* Convert the secs and usecs to floating point time */
  time = (float) secs + (float) usecs / 1000000.0F;
  avg_time = (float) avg_secs + (float) avg_usecs / 1000000.0F;

  /* do the output */
  snprintf(message, MAX_OUTSTR, "Calls: %" PRIdinteger, (integer)calls);
  umPrint(message,"pio_io_timer");

  snprintf(message, MAX_OUTSTR,
     "Total time: %.6f secs \t\t Average time: %.6f secs",
     time, avg_time );
  umPrint(message,"pio_io_timer");

  if (op == IO_B_IN || op == IO_B_OUT)
  {
    snprintf(message, MAX_OUTSTR,
           "Total bytes transferred: %lld \t Average bytes transferred:"
           " %lld", (long long int) bytes, (long long int) avg_bytes);
    umPrint(message,"pio_io_timer");

    /* Only print transfer rates if exists */
    if ( (time > 0.000001) && (avg_time > 0.000001) )
    {
      /* Note that transfers are specified in terms of Megabytes */
      /* per second in the SI sense - ie 10^6 bytes.             */
      /* This can easily be changed for Mebibytes = 2^20 bytes   */
      mbs = (float) bytes / time;
      mbs = mbs/ (1000.0F * 1000.0F);

      avg_mbs = (float) avg_bytes / avg_time;
      avg_mbs = avg_mbs/ (1000.0F * 1000.0F);

      snprintf(message, MAX_OUTSTR,
              "Overall transfer rate: %.6f MB/s \t Average transfer"
              " rate: %.6f MB/s", mbs, avg_mbs);
      umPrint(message,"pio_io_timer");
    }
  }

  return;
}

/****************************************************************/
/* The following function is to be called from Fortran.         */
/* It produces a simple report (on stdout) that details the     */
/* io usage for the complete program thus far.                  */
/****************************************************************/
void io_total_timings (void)
{
   char message[MAX_OUTSTR];

   /* BUFFIN */
   snprintf(message, MAX_OUTSTR, "\nTotal IO Timings: BUFFIN");
   umPrint(message,"pio_io_timer");

   io_report( io_calls[IO_B_IN], io_bytes[IO_B_IN], io_secs[IO_B_IN],
              io_usecs[IO_B_IN], IO_B_IN );


   /* BUFFOUT */
   snprintf(message, MAX_OUTSTR, "\nTotal IO Timings: BUFFOUT");
   umPrint(message,"pio_io_timer");

   io_report( io_calls[IO_B_OUT], io_bytes[IO_B_OUT], io_secs[IO_B_OUT],
              io_usecs[IO_B_OUT], IO_B_OUT );


   /* FILE_OPEN */
   snprintf(message, MAX_OUTSTR, "\nTotal IO Timings: FILE_OPEN");
   umPrint(message,"pio_io_timer");

   io_report( io_calls[IO_F_OP], io_bytes[IO_F_OP], io_secs[IO_F_OP],
              io_usecs[IO_F_OP], IO_F_OP );


   /* FILE_CLOSE */
   snprintf(message, MAX_OUTSTR, "\nTotal IO Timings: FILE_CLOSE");
   umPrint(message,"pio_io_timer");

   io_report( io_calls[IO_F_CL], io_bytes[IO_F_CL], io_secs[IO_F_CL],
              io_usecs[IO_F_CL], IO_F_CL );

   /* FILE_SEEK */
   snprintf(message, MAX_OUTSTR, "\nTotal IO Timings: FILE_SEEK");
   umPrint(message, "pio_io_timer");

   io_report(io_calls[IO_F_SK], io_bytes[IO_F_SK], io_secs[IO_F_SK],
             io_usecs[IO_F_SK], IO_F_SK);

   return;
}

/****************************************************************/
/* The following function is to be called from Fortran.         */
/* It produces a simple report (on stdout) that details the     */
/* io usage since the last call of this function.               */
/****************************************************************/
void io_output_inter_timings(void)
{
   char message[MAX_OUTSTR];
   int i;

   /* BUFFIN */
   snprintf(message, MAX_OUTSTR, "\nIntermediate IO Timings: BUFFIN");
   umPrint(message,"pio_io_timer");

   io_report( io_calls_last[IO_B_IN], io_bytes_last[IO_B_IN],
              io_secs_last[IO_B_IN], io_usecs_last[IO_B_IN], IO_B_IN );


   /* BUFFOUT */
   snprintf(message, MAX_OUTSTR, "\nIntermediate IO Timings: BUFFOUT");
   umPrint(message,"pio_io_timer");

   io_report( io_calls_last[IO_B_OUT], io_bytes_last[IO_B_OUT],
              io_secs_last[IO_B_OUT], io_usecs_last[IO_B_OUT],
              IO_B_OUT );


   /* FILE_OPEN */
   snprintf(message, MAX_OUTSTR, "\nIntermediate IO Timings: FILE_OPEN");
   umPrint(message,"pio_io_timer");

   io_report( io_calls_last[IO_F_OP], io_bytes_last[IO_F_OP],
              io_secs_last[IO_F_OP], io_usecs_last[IO_F_OP], IO_F_OP );


   /* FILE_CLOSE */
   snprintf(message, MAX_OUTSTR, "\nIntermediate IO Timings: FILE_CLOSE");
   umPrint(message,"pio_io_timer");

   io_report( io_calls_last[IO_F_CL], io_bytes_last[IO_F_CL],
              io_secs_last[IO_F_CL], io_usecs_last[IO_F_CL], IO_F_CL );

   /* reset the *last* data */
   for (i=1; i<4; i++)
   {
      io_calls_last[i] = 0;
      io_bytes_last[i] = 0;
      io_secs_last[i] = 0;
      io_usecs_last[i] = 0;
   }

   return;
}

/****************************************************************/
/* The following function is to be called from Fortran. It      */
/* is used to set the value of the io_timer_active flag.        */
/* 0        = disabled io_timer                                 */
/* non-zero = active io_timer                                   */
/****************************************************************/
void io_timing_init(void)
{
  io_timer_active = 1;
}

#endif
