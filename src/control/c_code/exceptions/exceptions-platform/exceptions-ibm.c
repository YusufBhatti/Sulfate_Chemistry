/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* Description:                                                               */
/*   C signal trapping routines,                                              */
/*   with a sensible default list of signals to trap.                         */

#if defined (IBM_XL_LIBC)

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif

#include <fptrap.h>
#include <stdio.h>
#include <inttypes.h>
#include "exceptions.h"
#include "exceptions-ibm.h"

/*----------------------------------------------------------------------------*/
/* Procedures                                                                 */
/*----------------------------------------------------------------------------*/

void signal_trap_extra_fp_ibm(fpetraps_t *set_traps)
{
  fptrap_t exceptions = 0;

  fp_trap(FP_TRAP_FASTMODE);

  #if defined(TRP_INVALID)
  if( set_traps->invalid )
  {
    exceptions |= TRP_INVALID;
  }
  #else
  if( set_traps->invalid )
  {
     fprintf(stderr,
             "Ignoring requested Floating-Point Invalid trapping"
             " as it is not supported on this platform");
  }
  #endif

  #if defined(TRP_DIV_BY_ZERO)
  if( set_traps->divbyzero )
  {
    exceptions |= TRP_DIV_BY_ZERO;
  }
  #else
  if( set_traps->divbyzero )
  {
     fprintf(stderr,
             "Ignoring requested Floating-Point Divide-by-Zero trapping"
             " as it is not supported on this platform");
  }
  #endif

  #if defined(TRP_OVERFLOW)
  if( set_traps->overflow )
  {
    exceptions |= TRP_OVERFLOW;
  }
  #else
  if( set_traps->overflow )
  {
     fprintf(stderr,
             "Ignoring requested Floating-Point Overflow trapping"
             " as it is not supported on this platform");
  }
  #endif

  #if defined(TRP_INEXACT)
  if( set_traps->inexact )
  {
    exceptions |= TRP_INEXACT;
  }
  #else
  if( set_traps->inexact )
  {
     fprintf(stderr,
             "Ignoring requested Floating-Point Inexact trapping"
             " as it is not supported on this platform");
  }
  #endif

  #if defined(TRP_UNDERFLOW)
  if( set_traps->underflow )
  {
    exceptions |= TRP_UNDERFLOW;
  }
  #else
  if( set_traps->underflow )
  {
     fprintf(stderr,
             "Ignoring requested Floating-Point Invalid trapping"
             " as it is not supported on this platform");
  }
  #endif

  fp_enable(exceptions);
}

/*----------------------------------------------------------------------------*/

/* IBM backtrace */
void signal_do_backtrace_ibm(int signo,
                             siginfo_t *sigcode,
                             void *sigcontextptr)
{

  /* on IBM the core dump also provides a traceback */
  if( sig_state.option == signal_core )
  {
    return;
  }

  if( sig_state.verbose )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions: *** XL_TRACE ***\n",
            sig_state.rank);
  }

  /* Gives a traceback and stops the program */
  xl__trce(signo, sigcode, sigcontextptr);

}

/*----------------------------------------------------------------------------*/

void signal_do_core_ibm(int signo,
                        siginfo_t *sigcode,
                        void *sigcontextptr)
{
  /* Give a traceback and dump core - then stop */

  if( sig_state.verbose )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions: *** XL_DUMP ***\n",
            sig_state.rank);
  }

  xl__trcedump(signo, sigcode, sigcontextptr);
}

/*----------------------------------------------------------------------------*/

#else
extern int unused_unit;
#endif
