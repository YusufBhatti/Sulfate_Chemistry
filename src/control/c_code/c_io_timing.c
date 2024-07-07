#if defined(C95_2B)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#include "c_io_layers.h"

#define CIO_THIS_LAYER CIO_TIMING

#include "c_io_nextlayer.h"

#include "c_io_internal.h"

#include "c_io_timing.h"

#include <float.h>

/* Functions in this file:
 *
 * Layer interfaces and functionality:
 *   init          : initialise timing structures and pass through
 *   fini          : report timing and pass through
 *   open          : time pass through call
 *   close         : time pass through call
 *   in            : time pass through call
 *   out           : time pass through call
 *   getpos        : time pass through call
 *   setpos        : time pass through call
 *   change_mode   : time pass through call
 *   sync          : time pass through call
 *   query         : time pass through call
 *
 * Implentation auxiliary routines:
 *   c_io_timer    : increments the time/counters - timer start/stop
 *   c_io_getTime  : provides the time.
 *
 */

/*----------------------------------------------------------------------------*/
/* Enumerators and Typedefs                                                   */
/*----------------------------------------------------------------------------*/

enum opTypes     { ioStart, ioStop, ioPause, ioUnPause, numOpTypes };
enum opGroups    { ioOpen, ioClose, ioIn, ioOut, ioSeek, ioMisc, ioExtra,
                   numOpGroups };

/* Record keeping for a timer module... */

typedef struct c_io_timer_t {
  double time;
  size_t bytes;
  size_t count;
} c_io_timer_t;

typedef struct c_io_timer_dat_t {
  c_io_timer_t * timers;
} c_io_timer_dat_t;

typedef struct c_io_timer_extra_dat_t {
  void *next;
  double time;
} c_io_timer_extra_dat_t;

/*----------------------------------------------------------------------------*/
/* Global Variables                                                           */
/*----------------------------------------------------------------------------*/

static c_io_timer_dat_t *c_io_timer_dat;
static bool c_io_timing_switch;
static c_io_timer_extra_dat_t **c_io_timer_extra_dat_list;
static c_io_timer_extra_dat_t **c_io_timer_extra_dat;

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

static void io_total_timings (void);

static void c_io_timer       (int64_t,
                              int64_t,
                              int64_t,
                              int64_t);

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_init (c_io_init_t *init_params){
  size_t i,j;

  c_io_init_next_layer(CIO_NEXT_LAYER);

  if (init_params->timing_switch == 1){
    c_io_timing_switch=true;
    c_io_printf(-1,"init: timing: Timing module is active");
  }
  else
    c_io_timing_switch=false;

  c_io_timer_dat =
                malloc(sizeof(*c_io_timer_dat)*(size_t)(init_params->max_unit));

  if (c_io_timing_switch)
  {
    c_io_timer_extra_dat_list = malloc((init_params->max_unit)
                                       * sizeof(*c_io_timer_extra_dat_list));

    c_io_timer_extra_dat = malloc((init_params->max_unit)
                                  * sizeof(*c_io_timer_extra_dat));

    for (i=(init_params->min_unit);i<(init_params->max_unit);i++)
    {
      c_io_timer_dat[i].timers =
                            malloc(sizeof(*c_io_timer_dat->timers)*numOpGroups);

      for (j=0;j<numOpGroups;j++)
      {
        c_io_timer_dat[i].timers[j].count=0;
        c_io_timer_dat[i].timers[j].time=0.0;
        c_io_timer_dat[i].timers[j].bytes=0;
      }

      c_io_timer_extra_dat_list[i] =
                                  malloc(sizeof(*c_io_timer_extra_dat_list[i]));
      c_io_timer_extra_dat[i]=c_io_timer_extra_dat_list[i];
      c_io_timer_extra_dat_list[i]->next=NULL;

    }
  }

  return C_IO_INIT(init_params);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_fini (void)
{
  int64_t res;

  res=C_IO_FINI();

  io_total_timings();

  if (c_io_timing_switch)
  {
    size_t i;
    c_io_timer_extra_dat_t *free_list;

#if defined(__INTEL_COMPILER)
    /* This loop cannot be vectorised due to the presence of free() calls.
     * Suppress any vectorisation attempt with a novector #pragma, to ensure no
     * compiler warnings are generated about failed vectorisation attempts.
     */
   #pragma novector
#endif
    for (i=(params->min_unit);i<(params->max_unit);i++)
    {
      free(c_io_timer_dat[i].timers);

      free_list=c_io_timer_extra_dat_list[i];

      while (free_list->next!=NULL)
      {
        c_io_timer_extra_dat[i]=free_list->next;
        free(free_list);
        free_list=c_io_timer_extra_dat[i];
      }

      free(free_list);
    }

    free(c_io_timer_extra_dat_list);
    free(c_io_timer_extra_dat);
  }

  free(c_io_timer_dat);

  return res;
}

/*----------------------------------------------------------------------------*/

void io_total_timings (void)
{

#define DASHSTR ":----------------------------------------------"              \
                "-----------------------------------------------"

#define GBFACTOR (1024.0 * 1024.0 * 1024.0)
#define MBFACTOR (1024.0 * 1024.0)

  size_t i,j;
  size_t t_bytes,p_bytes;
  double t_time,p_time;
  bool extra=false;

  if (c_io_timing_switch)
  {
    c_io_printf(-1, DASHSTR);
    c_io_printf(-1,": %11s %6s %18s %18s %13s %13s %7s",
                "<-- Open ->",
                "Close",
                "<----------- In ->",
                "<---------- Out ->",
                "<---- Seek ->",
                "<---- Misc ->",
                "Rate");
    c_io_printf(-1, ":"
                " %4s %6s"
                " %6s"
                " %5s %5s %6s"
                " %5s %5s %6s"
                " %6s %6s"
                " %6s %6s"
                " %7s",
                "CNT","T",
                "T",
                "CNT","GB","T",
                "CNT","GB","T",
                "CNT","T",
                "CNT","T",
                "MB/s");
    c_io_printf(-1, DASHSTR);

    p_time=0.0;
    p_bytes=0;
    for (i=(params->min_unit);i<(params->max_unit);i++)
    {
      if (c_io_timer_dat[i].timers[ioOpen].count>0)
      {
        t_time=0.0;
        t_bytes=0;
        for (j=0;j<numOpGroups;j++)
        {
          if (j==ioExtra)
          {
            continue;
          }

          t_time +=c_io_timer_dat[i].timers[j].time;
          t_bytes+=c_io_timer_dat[i].timers[j].bytes;
        }
        p_time+=t_time;
        p_bytes+=t_bytes;
        if (t_time>0.0)
        {
          c_io_printf((int64_t)i,
                      ": %4zu %6.2f %6.2f %5zu %5.2f %6.2f %5zu"
                      " %5.2f %6.2f %6zu %6.2f %6zu %6.2f %7.2f",
                      c_io_timer_dat[i].timers[ioOpen].count,
                      c_io_timer_dat[i].timers[ioOpen].time,
                      c_io_timer_dat[i].timers[ioClose].time,
                      c_io_timer_dat[i].timers[ioIn].count,
                      (double)(c_io_timer_dat[i].timers[ioIn].bytes)/GBFACTOR,
                      c_io_timer_dat[i].timers[ioIn].time,
                      c_io_timer_dat[i].timers[ioOut].count,
                      (double)(c_io_timer_dat[i].timers[ioOut].bytes)/GBFACTOR,
                      c_io_timer_dat[i].timers[ioOut].time,
                      c_io_timer_dat[i].timers[ioMisc].count,
                      c_io_timer_dat[i].timers[ioMisc].time,
                      c_io_timer_dat[i].timers[ioSeek].count,
                      c_io_timer_dat[i].timers[ioSeek].time,
                      ((double)(t_bytes)/MBFACTOR)/t_time
                      );
        }
        else
        {
          c_io_printf((int64_t)i,
                      ": %4zu %6.2f %6.2f %5zu %5.2f %6.2f %5zu"
                      " %5.2f %6.2f %6zu %6.2f %6zu %6.2f",
                      c_io_timer_dat[i].timers[ioOpen].count,
                      c_io_timer_dat[i].timers[ioOpen].time,
                      c_io_timer_dat[i].timers[ioClose].time,
                      c_io_timer_dat[i].timers[ioIn].count,
                      (double)(c_io_timer_dat[i].timers[ioIn].bytes)/GBFACTOR,
                      c_io_timer_dat[i].timers[ioIn].time,
                      c_io_timer_dat[i].timers[ioOut].count,
                      (double)(c_io_timer_dat[i].timers[ioOut].bytes)/GBFACTOR,
                      c_io_timer_dat[i].timers[ioOut].time,
                      c_io_timer_dat[i].timers[ioMisc].count,
                      c_io_timer_dat[i].timers[ioMisc].time,
                      c_io_timer_dat[i].timers[ioSeek].count,
                      c_io_timer_dat[i].timers[ioSeek].time
                      );
        }
      }
    }
    c_io_printf(-1, DASHSTR);
    if (p_time>0.0)
    {
      c_io_printf(-1,": Total data %.2fMB Total time %.2fs (Rate %.2fMB/s)",
                  (double)(p_bytes)/1024.0/1024.0,
                  p_time,
                  ((double)(p_bytes)/1024.0/1024.0)/p_time);
    }
    else
    {
      c_io_printf(-1,": Total data %.2fMB Total time %.2fs",
                  (double)(p_bytes)/1024.0/1024.0,
                  p_time);
    }

    c_io_printf(-1, DASHSTR);

    for (i=(params->min_unit);i<(params->max_unit);i++)
    {
      extra = (c_io_timer_dat[i].timers[ioExtra].count>0);
      if (extra)
      {
        break;
      }
    }

    if(extra)
    {
      c_io_timer_extra_dat_t *extra_dat_list;
      double max, min, ave;
      size_t *tot_count;

      c_io_printf(-1,": Extra Timing By Unit");
      c_io_printf(-1,": %7s %6s", "Count", "Time");

      tot_count = calloc(params->max_unit, sizeof(*tot_count));

      for (i=(params->min_unit);i<(params->max_unit);i++)
      {
        if (c_io_timer_dat[i].timers[ioExtra].count>0)
        {
          tot_count[i] = c_io_timer_dat[i].timers[ioExtra].count;

          c_io_printf((int64_t)i, ": %7" PRId64 " %6.2f",
                      c_io_timer_dat[i].timers[ioExtra].count,
                      c_io_timer_dat[i].timers[ioExtra].time
                      );
        }
      }

      c_io_printf(-1,":");
      c_io_printf(-1,": Extra Timings List");
      c_io_printf(-1,":");

      for (i=(params->min_unit);i<(params->max_unit);i++)
      {
        if (tot_count[i]>0)
        {
          extra_dat_list = c_io_timer_extra_dat_list[i];

          c_io_printf(-1,": Unit %" PRId64 "", i);
          c_io_printf(-1,": %7s %6s", "#", "Time");

          for (j=0;j<tot_count[i];j++)
          {
            double new_t = extra_dat_list->time;
            extra_dat_list = extra_dat_list->next;

            c_io_printf(-1, ": %7zu %6.2f",j+1,new_t);
          }

          c_io_printf(-1,":");

        }
      }

      for (i=(params->min_unit);i<(params->max_unit);i++)
      {
        if (tot_count[i]>0)
        {
          extra_dat_list = c_io_timer_extra_dat_list[i];

          ave = 0.0;
          max = -DBL_MAX;
          min = DBL_MAX;

          for (j=0;j<tot_count[i];j++)
          {
            double new_t = extra_dat_list->time;
            extra_dat_list = extra_dat_list->next;

            ave = ((double)j/(double)(j+1))*(ave-new_t) + new_t;
            min = (min < new_t) ? min : new_t;
            max = (max > new_t) ? max : new_t;

          }

          c_io_printf(-1,":");

          c_io_printf(-1,": Extra Timings Stats (Unit %" PRId64 ")", i);
          c_io_printf(-1,": %6s %6s %9s", "Min", "Max", "Mean");
          c_io_printf(-1,": %6.2f %6.2f %9.5f",min,max,ave);

          c_io_printf(-1,":");
        }
      }

      free(tot_count);

      c_io_printf(-1, DASHSTR);
    }

  }

}

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_open   (int64_t unit,
                            char *filename,
                            int64_t fileStatus,
                            int64_t fileMode)
{
  int64_t retVal;
  c_io_timer(unit,ioOpen,ioStart,0);
  retVal=C_IO_OPEN(unit,filename,fileStatus,fileMode);
  c_io_timer(unit,ioOpen,ioStop,0);
  return retVal;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_close  (int64_t unit)
{
  int64_t retVal;
  c_io_timer(unit,ioClose,ioStart,0);
  retVal=C_IO_CLOSE(unit);
  c_io_timer(unit,ioClose,ioStop,0);
  return retVal;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_in     (int64_t unit,
                            char *array,
                            int64_t len,
                            int64_t word_len)
{
  int64_t numRead;
  c_io_timer(unit,ioIn,ioStart,len*word_len);
  numRead=C_IO_IN(unit,array,len,word_len);
  c_io_timer(unit,ioIn,ioStop, numRead*word_len);
  return numRead;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_out    (int64_t unit,
                            char *array,
                            int64_t len,
                            int64_t word_len)
{
  int64_t numWritten;
  c_io_timer(unit,ioOut,ioStart,len*word_len);
  numWritten=C_IO_OUT(unit,array,len,word_len);
  c_io_timer(unit,ioOut,ioStop, numWritten*word_len);
  return  numWritten;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_getpos (int64_t unit)
{
  int64_t retVal;
  c_io_timer(unit,ioMisc,ioStart,0);
  retVal=C_IO_GETPOS(unit);
  c_io_timer(unit,ioMisc,ioStop,0);
  return retVal;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_setpos (int64_t unit,
                            int64_t pos,
                            int64_t word_len)
{
  int64_t retVal;
  c_io_timer(unit,ioSeek,ioStart,0);
  retVal=C_IO_SETPOS(unit,pos,word_len);
  c_io_timer(unit,ioSeek,ioStop,0);
  return retVal;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_change_mode(int64_t unit,int64_t newMode)
{
  int64_t retVal;
  c_io_timer(unit,ioMisc,ioStart,0);
  retVal=C_IO_CHANGE_MODE(unit,newMode);
  c_io_timer(unit,ioMisc,ioStop,0);
  return retVal;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_sync(int64_t unit)
{
  int64_t retVal;
  c_io_timer(unit,ioMisc,ioStart,0);
  retVal=C_IO_SYNC(unit);
  c_io_timer(unit,ioMisc,ioStop,0);
  return retVal;
}

/*----------------------------------------------------------------------------*/

int64_t c_io_timing_query(int64_t unit,
                          query_type_t query,
                          void *in_pack,
                          void *out_pack)
{
  int64_t retVal;
  c_io_timer(unit,ioMisc,ioStart,0);
  retVal=C_IO_QUERY(unit, query, in_pack, out_pack);
  c_io_timer(unit,ioMisc,ioStop,0);
  return retVal;
}

/*----------------------------------------------------------------------------*/

void c_io_timer(int64_t unit, int64_t group, int64_t op, int64_t ammount)
{
  static double tStart[numOpGroups];

  if (c_io_timing_switch)
  {
    if (op==ioStart)
    {
      tStart[group]=c_io_getTime();
    }
    else if (op==ioStop)
    {
      double add_time = c_io_getTime()-tStart[group];

      c_io_timer_dat[unit].timers[group].count++;
      c_io_timer_dat[unit].timers[group].time+=add_time;
      c_io_timer_dat[unit].timers[group].bytes+=(size_t)ammount;

      if (group==ioExtra)
      {
        c_io_timer_extra_dat[unit]->time=add_time;
        c_io_timer_extra_dat[unit]->next =
                                    malloc(sizeof(*c_io_timer_extra_dat[unit]));
        c_io_timer_extra_dat[unit]=c_io_timer_extra_dat[unit]->next;
        c_io_timer_extra_dat[unit]->next=NULL;
      }
    }
    else
    {
      c_io_printf(unit,"WARNING: c_io_timing: bad opCode");
    }
  }
}

/*----------------------------------------------------------------------------*/

void c_io_start_extra_timer(int64_t unit)
{
  c_io_timer(unit,ioExtra,ioStart,0);
}

/*----------------------------------------------------------------------------*/

void c_io_stop_extra_timer(int64_t unit)
{
  c_io_timer(unit,ioExtra,ioStop,0);
}

/*----------------------------------------------------------------------------*/

#endif
