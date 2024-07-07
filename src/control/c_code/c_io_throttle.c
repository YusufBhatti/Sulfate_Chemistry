#if defined(C95_2B) && defined(THROTTLE_IO)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* This unit is used to "throttle" system calls for I/O.
 *
 * It does this by creating a timer, which measures a minimum time between
 * successive calls. (This can be set using the IO_THROTTLE_SEC and 
 * IO_THROTTLE_NSEC environment variables.) If this ammount of time has not
 * elapsed at the point of a subsequent I/O call, the thread spools until it has 
 * before making that call.
 *
 * Currently, only read calls are throttled.
 *
 * Note: It is the call rate which is throttled, not the I/O bandwidth. i.e. a
 *       subsequent calls will be spooled regardless of the data size of its I/O
 *       transaction. 
 */

/* Compliance to Posix specification IEEE Std 1003.1-2001 or later is requred
 * by this unit in order to access some functionality (e.g. timer_settime() )
 */

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#include "c_io_layers.h"

#define CIO_THIS_LAYER CIO_THROTTLE

#include "c_io_nextlayer.h"

#include "c_io_internal.h"
#include "c_io_throttle.h"

#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <signal.h>
#include <time.h>

/*  Functions in this file:
 *
 *  Layer interfaces:
 *    init
 *    in
 *    open
 *    close
 *    setpos
 *    fini
 *    getpos
 *    out
 *    change_mode
 *    sync
 *    query
 *
 */

/*----------------------------------------------------------------------------*/

/* Clock to be used to limit call rate */
static timer_t throttle_clock;

/* Time structures required to set throttle_clock expiration times */
static struct itimerspec set_itrspec;
static const struct timespec null_tspec = {0,0};

/*----------------------------------------------------------------------------*/
/* Implentation routines                                                      */
/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_init(c_io_init_t *init_params)
{
  int timer_ret;
  const char *env_str=NULL;
  struct sigevent null_event;
  struct timespec set_tspec;

  c_io_init_next_layer(CIO_NEXT_LAYER);

  /* Stop a sigalarm raise (and therefore exit) on timer expirery */
  null_event.sigev_notify = SIGEV_NONE;

  /* but set up throttle timers */
  timer_ret = timer_create(CLOCK_MONOTONIC, &null_event, &throttle_clock);

  if (timer_ret != 0)
  {
    c_io_printf_stderr(-1, 
                       "init: c_io_throttle: failed to create throttle clock");
    return c_io_fail;
  }

  /* initialise the time intervals */
  set_tspec.tv_sec=0;
  set_tspec.tv_nsec=0;

  env_str=getenv("IO_THROTTLE_SEC");

  if (env_str != NULL)
  {
    time_t tmp_sec;
    tmp_sec = (time_t)strtol(env_str, NULL, 10);
    if (tmp_sec>0)
    {
      set_tspec.tv_sec=tmp_sec;
    }
  }

  env_str=getenv("IO_THROTTLE_NSEC");

  if (env_str != NULL)
  {
    long tmp_nsec;
    tmp_nsec = strtol(env_str, NULL, 10);
    if (tmp_nsec>0)
    {
      set_tspec.tv_nsec=tmp_nsec;
    }
  }

  if (set_tspec.tv_sec!=0 || set_tspec.tv_nsec!=0)
  {
    c_io_printf(-1, "init: c_io_throttle: Throttle time %" PRIdMAX ".%09ld sec",
                (intmax_t)set_tspec.tv_sec, set_tspec.tv_nsec);
  }

  set_itrspec.it_interval=null_tspec;
  set_itrspec.it_value=set_tspec;

  /* start the timer */
  timer_settime(throttle_clock, 0, &set_itrspec, NULL);

  return C_IO_INIT(init_params);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_in (int64_t unit,
                          char    *array,
                          int64_t len,
                          int64_t word_len)
{
  struct itimerspec ret_itrspec;

  timer_gettime(throttle_clock, &ret_itrspec);

  while (ret_itrspec.it_value.tv_sec > 0 || ret_itrspec.it_value.tv_nsec > 0)
  {
    /* We want to allow the OS scheduler to context switch to higher priority
     * threads, so need to sleep here. We want to re-schedule in the minimum 
     * ammount of time to keep as tightly as possible to the throttle_clock
     * (sleep time is rounded up; this is the time until the thread is 
     * re-scheduled, not necessarily the time until return, so the timer always
     * overruns slightly), so use the sleep(0) idiom (which re-schedules next
     * cycle).
     * nanosleep() is thread-safe; sleep() is not.
     * [usleep() is thread-safe, but deprecated by POSIX.] 
     */
    nanosleep(&null_tspec,NULL);
    timer_gettime(throttle_clock, &ret_itrspec);
  }

  /* reset the timer */
  timer_settime(throttle_clock, 0, &set_itrspec, NULL);

  return C_IO_IN(unit,array,len,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_setpos (int64_t unit,
                              int64_t pos,
                              int64_t word_len)
{
  struct itimerspec ret_itrspec;

  timer_gettime(throttle_clock, &ret_itrspec);

  while (ret_itrspec.it_value.tv_sec > 0 || ret_itrspec.it_value.tv_nsec > 0)
  {
    /* We want to allow the OS scheduler to context switch to higher priority
     * threads, so need to sleep here. We want to re-schedule in the minimum 
     * ammount of time to keep as tightly as possible to the throttle_clock
     * (sleep time is rounded up; this is the time until the thread is 
     * re-scheduled, not necessarily the time until return, so the timer always
     * overruns slightly), so use the sleep(0) idiom (which re-schedules next
     * cycle).
     * nanosleep() is thread-safe; sleep() is not.
     * [usleep() is thread-safe, but deprecated by POSIX.] 
     */
    nanosleep(&null_tspec,NULL);
    timer_gettime(throttle_clock, &ret_itrspec);
  }

  /* reset the timer */
  timer_settime(throttle_clock, 0, &set_itrspec, NULL);

  return C_IO_SETPOS(unit,pos,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_getpos (int64_t unit)
{
  struct itimerspec ret_itrspec;

  timer_gettime(throttle_clock, &ret_itrspec);

  while (ret_itrspec.it_value.tv_sec > 0 || ret_itrspec.it_value.tv_nsec > 0)
  {
    /* We want to allow the OS scheduler to context switch to higher priority
     * threads, so need to sleep here. We want to re-schedule in the minimum 
     * ammount of time to keep as tightly as possible to the throttle_clock
     * (sleep time is rounded up; this is the time until the thread is 
     * re-scheduled, not necessarily the time until return, so the timer always
     * overruns slightly), so use the sleep(0) idiom (which re-schedules next
     * cycle).
     * nanosleep() is thread-safe; sleep() is not.
     * [usleep() is thread-safe, but deprecated by POSIX.] 
     */
    nanosleep(&null_tspec,NULL);
    timer_gettime(throttle_clock, &ret_itrspec);
  }

  /* reset the timer */
  timer_settime(throttle_clock, 0, &set_itrspec, NULL);

  return C_IO_GETPOS(unit);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_fini(void)
{
  return C_IO_FINI();
}

/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_change_mode (int64_t unit,
                                   int64_t newMode)
{
  return C_IO_CHANGE_MODE(unit,newMode);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_query (int64_t      unit,
                             query_type_t query,
                             void         *in_pack,
                             void         *out_pack)
{
  return C_IO_QUERY(unit, query, in_pack, out_pack);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_close (int64_t unit)
{
  return C_IO_CLOSE(unit);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_open (int64_t unit,
                            char    *filename,
                            int64_t fileStatus,
                            int64_t fileMode)
{
  return C_IO_OPEN(unit,filename,fileStatus,fileMode);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_out (int64_t unit,
                           char *array,
                           int64_t len,
                           int64_t word_len)
{
  return C_IO_OUT(unit,array,len,word_len);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_throttle_sync(int64_t unit)
{
  return C_IO_SYNC(unit);
}

/*----------------------------------------------------------------------------*/

#endif
