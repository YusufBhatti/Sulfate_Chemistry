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

#if !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif

#if !defined(_BSD_SOURCE)
#define _BSD_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "exceptions.h"
#include "exceptions-generic.h"

void signal_do_backtrace_generic(int signo,
                                 siginfo_t *sigcode,
                                 void *sigcontextptr)
{
  /* The arguments are not used by this platform variant, kept for
   * compatibility with other variants.
   *
   * We therefore void cast those variables here to avoid failures of
   * compiler "unused variable" checks: -
   */
  (void)(signo);
  (void)(sigcode);
  (void)(sigcontextptr);

  fprintf(stderr,
          "[%" PRId64 "] exceptions:"
          " No backtrace function defined for this platform\n",
          sig_state.rank);


}

void signal_do_core_generic(int signo,
                            siginfo_t *sigcode,
                            void *sigcontextptr)
{
  /* The arguments are not used by this platform variant, kept for
   * compatibility with other variants.
   *
   * We therefore void cast those variables here to avoid failures of
   * compiler "unused variable" checks: -
   */
  (void)(signo);
  (void)(sigcode);
  (void)(sigcontextptr);

  if (sig_state.verbose)
  {
      fprintf(stderr,
              "[%" PRId64 "] exceptions: *** ABORT ***\n",
              sig_state.rank);
  }

  abort();
}
