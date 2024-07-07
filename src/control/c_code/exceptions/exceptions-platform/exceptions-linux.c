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

#if defined (GNU_LIBC)

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif

#if !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include <fenv.h>
#include <stdio.h>
#include <signal.h>
#include <execinfo.h>
#include <ucontext.h>
#include <inttypes.h>
#include <stdint.h>
#include "exceptions.h"
#include "exceptions-linux.h"

/*----------------------------------------------------------------------------*/
/* Procedures                                                                 */
/*----------------------------------------------------------------------------*/

/* GLIBC backtrace */
void signal_do_backtrace_linux(int signo,
                               siginfo_t *sigcode,
                               void *sigcontextptr)
{
  void *array[1+MAX_TRACEBACK];
  int size;
  int i,j;
  void *addr=NULL;

  /* for file line info */
  char file[MAX_FILENAME];
  char func[MAX_FILENAME];
  int line;

  typedef struct sigcontext sigctxt_t;

  /* we don't need "signo" for this platform,
   * so avoid failures of compiler "unused variable" checks by void casting
   */
  (void)signo;

  ucontext_t *uc = (ucontext_t *)sigcontextptr;
  sigctxt_t  *sc = (sigctxt_t *)&(uc->uc_mcontext);

  addr = (void *)sc->INSTR_PTR;

  if( sig_state.verbose )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions: *** GLIBC ***\n",
            sig_state.rank);
  }

  fprintf(stderr,
          "[%" PRId64 "] exceptions:"
          " Data address (si_addr): 0x%08" PRIxPTR "; %s: 0x%08" PRIxPTR "\n",
          sig_state.rank, (uintptr_t)sigcode->si_addr,
          INSTR_PTR_STR, (uintptr_t)addr);

  size = backtrace(&array[1], MAX_TRACEBACK);
  array[0] = addr;

  fprintf(stderr,
          "[%" PRId64 "]"
          " exceptions: [backtrace]: has %3d elements:\n",
          sig_state.rank, size+1);

  for(i = 0; i <= size; i++)
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " [backtrace]: (%3d) : Address: [0x%08" PRIxPTR "] \n",
            sig_state.rank, i+1, (uintptr_t)array[i]);

    j = lineInfo(array[i], file, MAX_FILENAME, func, MAX_FILENAME, &line);

    switch(j)
    {
      case(0):
      {
        fprintf(stderr,
                "[%" PRId64 "] exceptions:"
                " [backtrace]: (%3d) : %s in file %s line %d\n",
                sig_state.rank, i+1, func, file, line);
      }
      break;

      case(2):
      {
        fprintf(stderr,
                "[%" PRId64 "] exceptions:"
                " [backtrace]: (%3d) : %s (%s)\n",
                sig_state.rank, i+1, func, file);
      }
      break;

      default:
      {
        fprintf(stderr,
                "[%" PRId64 "] exceptions: [backtrace]: (%3d) :"
                " line information unavailable error code: %d (%s)\n",
                sig_state.rank, i+1, j, file);
      }
      break;

    }

  }

  fprintf(stderr,
          "[%" PRId64 "] exceptions: \n",
          sig_state.rank);
  fprintf(stderr,
          "[%" PRId64 "] exceptions: "
          "To find the source line for an entry in the backtrace;\n",
          sig_state.rank);
  fprintf(stderr,
          "[%" PRId64 "] exceptions: "
          "run addr2line --exe=</path/too/executable> <address>\n",
          sig_state.rank);
  fprintf(stderr,
          "[%" PRId64 "] exceptions: "
          "where address is given as [0x<address>] above\n",
          sig_state.rank);
  fprintf(stderr,
          "[%" PRId64 "] exceptions: \n",
          sig_state.rank);

}

/*----------------------------------------------------------------------------*/

void signal_trap_extra_fp_linux(fpetraps_t *set_traps)
{
  int exceptions = 0;

  fedisableexcept(FE_ALL_EXCEPT);

  #if defined(FE_INVALID)
  if( set_traps->invalid )
  {
    exceptions |= FE_INVALID;
  }
  #else
  if( set_traps->invalid )
  {
     fprintf(stderr,
             "Ignoring requested Floating-Point Invalid trapping"
             " as it is not supported on this platform");
  }
  #endif

  #if defined(FE_DIVBYZERO)
  if( set_traps->divbyzero )
  {
    exceptions |= FE_DIVBYZERO;
  }
  #else
  if( set_traps->divbyzero )
  {
     fprintf(stderr,
             "Ignoring requested Floating-Point Divide-by-Zero trapping"
             " as it is not supported on this platform");
  }
  #endif

  #if defined(FE_OVERFLOW)
  if( set_traps->overflow )
  {
    exceptions |= FE_OVERFLOW;
  }
  #else
  if( set_traps->overflow )
  {
     fprintf(stderr,
             "Ignoring requested Floating-Point Overflow trapping"
             " as it is not supported on this platform");
  }
  #endif

  #if defined(FE_INEXACT)
  if( set_traps->inexact )
  {
    exceptions |= FE_INEXACT;
  }
  #else
  if( set_traps->inexact )
  {
     fprintf(stderr,
             "Ignoring requested Floating-Point Inexact trapping"
             " as it is not supported on this platform");
  }
  #endif

  #if defined(FE_UNDERFLOW)
  if( set_traps->underflow )
  {
    exceptions |= FE_UNDERFLOW;
  }
  #else
  if( set_traps->underflow )
  {
    fprintf(stderr,
            "Ignoring requested Floating-Point Invalid trapping"
            " as it is not supported on this platform");
  }
  #endif

  feenableexcept(exceptions);
  fprintf(stderr, "[%" PRId64 "] exceptions:"
          " feenableexcept() mask [0x%08x] enabled."
          " (mask [0x%08x] requested)\n",
          sig_state.rank, fegetexcept(), exceptions);

}

/*----------------------------------------------------------------------------*/

#else
extern int unused_unit;
#endif
