#if defined(__linux__)
/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* C language routines used to determine affinity information.                */

#if !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include <sched.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/syscall.h>

/*------------------------------------------------------------------------------
* Prototypes
*-------------------------------------------------------------------------------
*/

int c_max_available_cpus(void);
int c_num_available_cpus(void);
int c_running_on_core(void);
int cpumask_weight(cpu_set_t *cpumask);

/*------------------------------------------------------------------------------
* SYNOPSIS
*   int c_max_available_cpus()
*
* DESCRIPTION
*   Returns the maximum number of cores (real and virtual) available on a node.
*-------------------------------------------------------------------------------
*/

int c_max_available_cpus(void)
{
  int max_cpus=999;
  max_cpus = (int) sysconf(_SC_NPROCESSORS_ONLN);
  return max_cpus;
}

/*------------------------------------------------------------------------------
* SYNOPSIS
*   int num_available_cpus()
*
* DESCRIPTION
*   Returns the number of cores (real and virtual) on which this thread may run.
*-------------------------------------------------------------------------------
*/

int c_num_available_cpus(void)
{
  int num_cpus;
  cpu_set_t cpumask;
  pid_t tid;

  tid =  (pid_t) syscall(SYS_gettid);
  CPU_ZERO(&cpumask);

  sched_getaffinity(tid, sizeof(cpu_set_t), &cpumask);

  num_cpus = cpumask_weight(&cpumask);

  return num_cpus;
}

/*------------------------------------------------------------------------------
* SYNOPSIS
*   int cpumask_weight()
*
* DESCRIPTION
*   Returns the number of elements of the cpumask that are set.
*-------------------------------------------------------------------------------
*/

int cpumask_weight(cpu_set_t * cpumask)
{
  int i;
  int w;

  w=0;
  for (i=0; i < CPU_SETSIZE; i++){
    if (CPU_ISSET( (size_t) i , cpumask)){w++;}
  }

  return w;
}

/*------------------------------------------------------------------------------
* SYNOPSIS
*   int c_running_on_core()
*
* DESCRIPTION
*   Returns the ID of the core on which the calling task/thread is running.
*-------------------------------------------------------------------------------
*/

int c_running_on_core(void)
{
  int core=999;
  core = sched_getcpu();
  return core;
}

#endif

