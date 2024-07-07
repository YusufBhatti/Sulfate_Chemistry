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

#if defined (LIBUNWIND)

#include <libunwind.h>
#include <stdio.h>
#include <inttypes.h>
#include "exceptions.h"
#include "exceptions-libunwind.h"

/*----------------------------------------------------------------------------*/
/* Procedures                                                                 */
/*----------------------------------------------------------------------------*/

/* libunwind backtrace */
void signal_do_backtrace_libunwind(int signo,
                                   siginfo_t *sigcode,
                                   void *sigcontextptr)
{
  unw_cursor_t cursor;
  unw_context_t uc;
  unw_word_t ip, sp;
  char bufp[MAX_FILENAME];
  unw_word_t off;
  int i,j,error;

  /* for file line info */
  char file[MAX_FILENAME];
  char func[MAX_FILENAME];
  int line;

  if( sig_state.verbose )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions: *** LIBUNWIND ***\n",
            sig_state.rank);
  }

  unw_getcontext(&uc);
  unw_init_local(&cursor, &uc);

  if( sig_state.verbose )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions: *** LIBUNWIND TRACEBACK ***\n",
            sig_state.rank);
  }

  i=1;

  while( (error=unw_step(&cursor)) > 0 )
  {
    if( unw_get_reg(&cursor, UNW_REG_IP, &ip) != 0 )
    {
      fprintf(stderr,
          "[%" PRId64 "] exceptions: [backtrace]: (%3d) : "
          "Recovering the instruction pointer failed\n",
          sig_state.rank, i);
    }

    if( unw_get_reg(&cursor, UNW_REG_SP, &sp) != 0 )
    {
      fprintf(stderr,
              "[%" PRId64 "] exceptions: [backtrace]: (%3d) : "
              "Recovering the stack pointer failed\n",
              sig_state.rank, i);
    }

    if( unw_is_signal_frame(&cursor) > 0 )
    {
      fprintf(stderr,
              "[%" PRId64 "] exceptions: [backtrace]: (%3d) : "
              "*v* SIGNAL GENERATING FRAME *v*\n",
              sig_state.rank, i);
    }

    if( (j=unw_get_proc_name(&cursor,bufp,MAX_FILENAME,&off)) != 0 )
    {
      if( j == UNW_EUNSPEC )
      {
        fprintf(stderr,
                "[%" PRId64 "] exceptions: [backtrace]: (%3d) :  "
                "Obtaining the procedure name/byte offset failed "
                "(UMW_EUNSPEC)\n",
                sig_state.rank, i);
      }
      else if( j == UNW_ENOINFO )
      {
        fprintf(stderr,
                "[%" PRId64 "] exceptions: [backtrace]: (%3d) :  "
                "Obtaining the procedure name/byte offset failed "
                "(UMW_ENOINFO)\n",
                sig_state.rank, i);
      }
      else if( j == UNW_ENOINFO )
      {
        fprintf(stderr,
                "[%" PRId64 "] exceptions: [backtrace]: (%3d) :  "
                "Obtaining the procedure name/byte offset failed "
                "(UNW_ENOMME)\n",
                sig_state.rank, i);
      }
      else
      {
        fprintf(stderr,
                "[%" PRId64 "] exceptions: [backtrace]: (%3d) :  "
                "Obtaining the procedure name/byte offset failed "
                "(Unknown %d)\n",
                sig_state.rank, i, j);
      }
    }

    fprintf(stderr,
            "[%" PRId64 "] exceptions: [backtrace]: (%3d) : "
            "%s + %" PRIuMAX " ip=0x%08" PRIxMAX ", sp=0x%08" PRIxMAX " \n",
            sig_state.rank, i, bufp, (uintmax_t)off,
            (uintmax_t)ip, (uintmax_t)sp);

    j = lineInfo(ip, file, MAX_FILENAME, func, MAX_FILENAME, &line);

    if( j == 0 )
    {
      fprintf(stderr,
              "[%" PRId64 "] exceptions: [backtrace]: (%3d) : "
              "%s in file %s line %d\n",
              sig_state.rank, i, func, file, line);
    }
    else
    {
      fprintf(stderr,
              "[%" PRId64 "] exceptions: [backtrace]: (%3d) : "
              "line information unavailable error code: %d (%s)\n",
              sig_state.rank, i, j, file);
    }

    if( unw_is_signal_frame(&cursor) > 0 )
    {
      fprintf(stderr,
              "[%" PRId64 "] exceptions: [backtrace]: (%3d) : "
              "*^* SIGNAL GENERATING FRAME *^*\n",
              sig_state.rank, i);
    }

    i++;
  }

  if( error == 0 )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " *** UNWINDING THE STACK ENDED OK ***\n",
            sig_state.rank);
  }
  else if( error < 0 )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " *** UNWINDING THE STACK ENDED BADLY ***\n",
            sig_state.rank);
  }

  if( sig_state.verbose )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " *** LIBUNWIND TRACEBACK DONE ***\n",
            sig_state.rank);
  }
}

/*----------------------------------------------------------------------------*/

#else
extern int unused_unit;
#endif
