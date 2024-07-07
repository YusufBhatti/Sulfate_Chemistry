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

#include <fenv.h>
#include <stdint.h>
#include <memory.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <inttypes.h>
#include <signal.h>
#include <stdint.h>

#include "thread_utils.h"
#include "exceptions.h"
#include "exceptions-generic.h"

/* Platform specific includes */

#if defined(GNU_LIBC)
#include "exceptions-linux.h"
#endif

#if defined(LIBUNWIND)
#include "exceptions-libunwind.h"
#endif

#if defined(IBM_XL_LIBC)
#include "exceptions-ibm.h"
#endif

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

void signal_trap                   (int64_t);
void signal_trap_extra             (void);
void signal_handler                (int, siginfo_t *, void *);
void signal_rank                   (int64_t pe);
void signal_set_verbose            (void);
void signal_command                (char *, int64_t);
void signal_do_core                (int, siginfo_t *, void *);
void signal_do_backtrace           (int, siginfo_t *, void *);
void signal_do_callbacks           (void);
void signal_controlled_exit        (int64_t);
void signal_add_handler            (void (*func)(void));
void signal_report_circumstance    (void);
void signal_report_exception_extra (int, siginfo_t *, int64_t);
void signal_unregister_callbacks   (void);

/*----------------------------------------------------------------------------*/
/* Global variables                                                           */
/*----------------------------------------------------------------------------*/

/* These are the unix signals we will intercept when active */
static int default_sigs[] = {
                             #if defined(SIGTRAP)
                             SIGTRAP,
                             #endif
                             #if defined(SIGFPE)
                             SIGFPE,
                             #endif
                             #if defined(SIGILL)
                             SIGILL,
                             #endif
                             #if defined(SIGBUS)
                             SIGBUS,
                             #endif
                             #if defined(SIGSEGV)
                             SIGSEGV,
                             #endif
                             #if defined(SIGXCPU)
                             SIGXCPU,
                             #endif
                             0
                            };

/* allows only one traceback to occur */
static volatile sig_atomic_t signal_trap_incr = 0;


/* global instance of state with initialisation */
signal_state_t sig_state = {             /* default to: -                     */
                            signal_off,  /*   don't handle signals            */
                            0,           /*   quiet output                    */
                            NULL,        /*   no callbacks                    */
                            0,           /*   unititialised processor id      */
                            0,           /*   unititialised MPI rank          */
                            "",          /*   unititialised program name      */
                            ""           /*   unititialised hostname          */
                           };

static bool extra_fp_options = false;

static fpetraps_t set_fpe_traps = {        /* default to: -                   */
                                   false,  /*   don't trap division by zero   */
                                   false,  /*   don't trap FP invalid         */
                                   false,  /*   don't trap overflow           */
                                   false,  /*   don't trap underflow          */
                                   false   /*   don't trap inexact results    */
                                  };

/*----------------------------------------------------------------------------*/
/* Procedures                                                                 */
/*----------------------------------------------------------------------------*/

/*
 * Generate a corefile if required.
 * This most likely stops the program too
 */
void signal_do_core(int signo,
                    siginfo_t *sigcode,
                    void *sigcontextptr)
{
  if (sig_state.option == signal_core)
  {
    #if defined (IBM_XL_LIBC)

    signal_do_core_ibm(signo, sigcode, sigcontextptr);

    #else

    signal_do_core_generic(signo, sigcode, sigcontextptr);

    #endif
  }
}

/*----------------------------------------------------------------------------*/

/* recover file/line info from the instruction pointer */
int lineInfo(void *addr,
             char *file,
             size_t filelen,
             char *func,
             size_t funclen,
             int *line)
{
  static char buf[MAX_FILENAME];
  char *p;
  int errCode=0; /* sucess */
  char *newline_ptr=NULL;
  size_t basewidth=0;

#define LINEINFO_ADDR2LINE_CMD "addr2line -C -e %.*s -f -i %lx"

  if( strlen(sig_state.command) < 1 )
  {
    snprintf(file, filelen, "Program name unknown");
    errCode=3;
    return errCode;
  }

  if( sig_state.command[0] != '/' )
  {
    snprintf(file, filelen, "Program name & path is not absolute");
    errCode=4;
    return errCode;
  }

  snprintf(buf, MAX_FILENAME, LINEINFO_ADDR2LINE_CMD,
           0, "", (unsigned long)addr);

  basewidth = strlen(buf);

  if( (basewidth + strlen(sig_state.command))  > MAX_FILENAME-1 )
  {
    errCode=5;
    snprintf(file, filelen, "Program name & path are too long");
    return errCode;
  }

  snprintf(buf, MAX_FILENAME, LINEINFO_ADDR2LINE_CMD,
           (int)(MAX_FILENAME-basewidth-1), sig_state.command,
           (unsigned long)addr);

  FILE *f = popen(buf, "r");

  if(f == NULL)
  {
    perror (buf);
    errCode=1;
    snprintf(file,filelen,"popen failed");
    return errCode;
  }

  /* find the function name */
  if( fgets(buf, MAX_FILENAME, f) == NULL )
  {
    snprintf(file, filelen, "Unable to get function name");
    return 10;
  }
  strncpy(func, buf, funclen-1);

  newline_ptr = strchr(func, '\n');

  if( newline_ptr != NULL )
  {
    *newline_ptr = '\0';
  }

  /* find the filename */
  if( fgets(buf, MAX_FILENAME, f) == NULL)
  {
    snprintf(file, filelen, "Unable to get file name");
    return 20;
  }

  if( buf[0] != '?' )
  {
    /* put a \0 in buff over the 1st colon */
    p = buf;

    while(*p != ':')
    {
      p++;
    }

    *p++ = 0;

    /* and rescue the filename and line number */
    strncpy(file, buf, filelen-1);
    sscanf(p, "%d", line);
  }
  else
  {
    /* fail sensibly */
    strncpy (file,"* Cannot Locate *", filelen-1);
    *line = 0;
    errCode=2;
  }

  pclose(f);

  return errCode;
}

/*----------------------------------------------------------------------------*/

/* backtrace - switch between compiled in choices */
void signal_do_backtrace(int signo,
                         siginfo_t *sigcode,
                         void *sigcontextptr)
{
  if( sig_state.option != signal_off )
  {
    #if defined(LIBUNWIND)

    signal_do_backtrace_libunwind(signo,sigcode,sigcontextptr);

    #elif defined(GNU_LIBC)

    signal_do_backtrace_linux(signo,sigcode,sigcontextptr);

    #elif defined(IBM_XL_LIBC)

    signal_do_backtrace_ibm(signo,sigcode,sigcontextptr);

    #else

    signal_do_backtrace_generic(signo,sigcode,sigcontextptr);

    #endif

  }
}

/*----------------------------------------------------------------------------*/

/* Call all of the registed callbacks */
void signal_do_callbacks(void)
{
  signal_callback_t *next;

  if( sig_state.callbacks == NULL )
  {
    if( sig_state.verbose )
    {
      fprintf(stderr,
              "[%" PRId64 "] exceptions: No callbacks \n",
              sig_state.rank);
    }
  }
  else
  {
    next=sig_state.callbacks;

    while( next != NULL )
    {
      if( sig_state.verbose )
      {
        fprintf(stderr,
                "[%" PRId64 "] exceptions:"
                " calling registered handler @ 0x%08" PRIxPTR "\n",
                sig_state.rank,
                (uintptr_t)next->callback);
      }

      next->callback();
      next=next->next;
    }
  }
}

/*----------------------------------------------------------------------------*/

/* Unregisted all of the callbacks */
void signal_unregister_callbacks(void)
{
  signal_callback_t *next=sig_state.callbacks;
  signal_callback_t *tmp;

  while (next != NULL)
  {
      tmp=next->next;
      free(next);
      next=tmp;
  }
}

/*----------------------------------------------------------------------------*/

/* What we know: extra information */
void signal_report_exception_extra(int signo,
                                   siginfo_t *sigcode,
                                   int64_t rank)
{
  const char *extrastr = NULL;

  /* deal with signal specific information */
  switch(signo)
  {
    #if defined(SIGFPE)
    case SIGFPE:
    {
      switch(sigcode->si_code)
      {
        #if defined(FPE_INTDIV)
        case FPE_INTDIV:
        {
          extrastr = "Integer divide by zero.";
        }
        break;
        #endif

        #if defined(FPE_INTOVF)
        case FPE_INTOVF:
        {
          extrastr = "Integer overflow.";
        }
        break;
        #endif

        #if defined(FPE_FLTDIV)
        case FPE_FLTDIV:
        {
          extrastr = "Floating-point divide by zero.";
        }
        break;
        #endif

        #if defined(FPE_FLTOVF)
        case FPE_FLTOVF:
        {
          extrastr = "Floating-point overflow.";
        }
        break;
        #endif

        #if defined(FPE_FLTUND)
        case FPE_FLTUND:
        {
          extrastr = "Floating-point underflow.";
        }
        break;
        #endif

        #if defined(FPE_FLTRES)
        case FPE_FLTRES:
        {
          extrastr = "Floating-point inexact result.";
        }
        break;
        #endif

        #if defined(FPE_FLTINV)
        case FPE_FLTINV:
        {
          extrastr = "Floating-point invalid operation.";
        }
        break;
        #endif

        #if defined(FPE_FLTSUB)
        case FPE_FLTSUB:
        {
          extrastr = "Subscript out of range.";
        }
        break;
        #endif

        default:
        break;
      }
    }
    break;
    #endif

    #if defined(SIGSEGV)
    case SIGSEGV:
    {
      switch(sigcode->si_code)
      {
        #if defined(SEGV_MAPERR)
        case SEGV_MAPERR:
        {
          extrastr = "Address not mapped to object.";
        }
        break;
        #endif

        #if defined(SEGV_ACCERR)
        case SEGV_ACCERR:
        {
          extrastr = "Invalid permissions for mapped object.";
        }
        break;
        #endif

        default:
        break;
      }
    }
    #endif

    default:
    break;
  }

  /* deal with generic signal codes */
  switch(sigcode->si_code)
  {
    #if defined(SI_USER)
    case SI_USER:
    {
      extrastr = "Sent by kill(2).";
    }
    break;
    #endif

    #if defined(SI_KERNEL)
    case SI_KERNEL:
    {
      extrastr = "Sent by the kernel.";
    }
    break;
    #endif

    #if defined(SI_QUEUE)
    case SI_QUEUE:
    {
      extrastr = "Sent by sigqueue(3).";
    }
    break;
    #endif

    #if defined(SI_TIMER)
    case SI_TIMER:
    {
      extrastr = "POSIX timer expired.";
    }
    break;
    #endif

    #if defined(SI_MESGQ)
    case SI_MESGQ:
    {
      extrastr = "POSIX message queue state changed";
    }
    break;
    #endif

    #if defined(SI_ASYNCIO)
    case SI_ASYNCIO:
    {
      extrastr = "AIO completed.";
    }
    break;
    #endif

    #if defined(SI_SIGIO)
    case SI_SIGIO:
    {
      extrastr = "Queued SIGIO";
    }
    break;
    #endif

    #if defined(SI_TKILL)
    case SI_TKILL:
    {
      extrastr = "Sent by tkill(2) or tgkill(2).";
    }
    break;
    #endif

    default:
    break;
  }

  /* If we found some extra information, print it now */
  if( extrastr != NULL )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " the exception reports the extra information: %s\n",
            rank, extrastr);
  }

}

/*----------------------------------------------------------------------------*/

/* What we know */
void signal_report_circumstance(void)
{
  if( inPar() == 1 )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " whilst in a parallel region, by thread %" PRId64 "\n",
            sig_state.rank,
            threadID());
  }
  else
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " whilst in a serial region\n",
            sig_state.rank);
  }

  fprintf(stderr,
          "[%" PRId64 "] exceptions:"
          " Task had pid=%" PRIdMAX " on host %s\n",
          sig_state.rank,
          (intmax_t)sig_state.pid,
          sig_state.hostname);
  fprintf(stderr,
          "[%" PRId64 "] exceptions: Program is \"%s\"\n",
          sig_state.rank,
          sig_state.command);
}

/*----------------------------------------------------------------------------*/

/*
 * for situations where we don't/can't call ereport
 *  we will do cleanup and produce traceback
 */
void signal_controlled_exit(int64_t e)
{
  fprintf(stderr,
          "[%" PRId64 "] exceptions:"
          " An non-exception application exit occured.\n",
          sig_state.rank);

  signal_report_circumstance();
  signal_do_callbacks();

  if( sig_state.verbose )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " Done callbacks\n",
            sig_state.rank);
  }

  signal_unregister_callbacks();

  if( e == 1 )
  {
    exit(255);
  }
}

/*----------------------------------------------------------------------------*/

/* The signal handler for fatal exception events */
void signal_handler(int signo,
                    siginfo_t *sigcode,
                    void *sigcontextptr)
{
  /*
   * We set an alarm which will terminate the app after the specified time.
   * Subsequent threads calling the handler will be punted into the long
   * grass for longer than that, so we don't need to worry about them.
   */

  /* Register the alarm to trigger after ALARM_CALL seconds  */
  alarm(ALARM_CALL);
  /* handle it with the system default handler */
  signal(SIGALRM, SIG_DFL);
  /* note: SIG_DFL is the default system handler, which is to terminate */

  /* Calling asynchronous interrupt unsafe functions from a signal handler
   * may lead to undefined behaviour or deadlock.
   *
   * From this point we will call such functions.
   *
   * In the event of a deadlock, we will be released by the Alarm Call.
   * The output following this point may however become corrupted.
   */

  if( signal_trap_incr++ )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " Subsequent thread (%" PRId64 ") calling handler will sleep\n",
            sig_state.rank,
            threadID());
    sleep(ALARM_CALL+10);
  }

  fprintf(stderr,
          "[%" PRId64 "] exceptions:"
          " An exception was raised:%d (%s)\n",
          sig_state.rank,
          signo,strsignal(signo));

  signal_report_exception_extra(signo, sigcode, sig_state.rank);
  signal_report_circumstance();

  if( sig_state.option == signal_core )
  {
    /* In systems with Large File Support, core_limit should be the 64-bit
     * type rlimit64, rather than the 32-bit rlimit.
     *
     * As an extension, some compilers (e.g. GCC with glibc) map rlimit to
     * rlimit64 automatically (and transparently) when
     * _FILE_OFFSET_BITS == 64.
     *
     * The Intel Compilers do not do this automatically.
     */
#if (_FILE_OFFSET_BITS == 64) && (__INTEL_COMPILER)
    struct rlimit64 core_limit;

    getrlimit64(RLIMIT_CORE, &core_limit);
#else
    struct rlimit core_limit;

    getrlimit(RLIMIT_CORE, &core_limit);
#endif

    if( core_limit.rlim_cur == 0 )
    {
      fprintf(stderr,
              "[%" PRId64 "] exceptions:"
              " Cannot core:  set \"ulimit -c\" to enable\n",
              sig_state.rank);
      sig_state.option = signal_traceback;
    }
  }

  signal_do_callbacks();

  /* it is unsafe to call signal_unregister_callbacks() here as it contains a
   * free() call. We will cause a memory leak [the memory malloc()'d to
   * sig_state.callbacks will never be free()'d], but that is preferable in
   * this instance to deadlocking or corrupting output.
   */

  if( sig_state.verbose )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " Done callbacks\n",
            sig_state.rank);
  }

  signal_do_backtrace(signo,sigcode,sigcontextptr);

  if( sig_state.verbose && (sig_state.option != signal_off) )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " Done backtrace\n",
            sig_state.rank);
  }

  signal_do_core(signo,sigcode,sigcontextptr);

  if( sig_state.verbose && (sig_state.option == signal_core) )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " Done core\n",
            sig_state.rank);
  }

  /* _exit() with the signo in situations where we reached the end */
  _exit(signo);
}

/*----------------------------------------------------------------------------*/

/* Set up exception handling */
void signal_trap(int64_t option)
{
  struct sigaction sa;

  int  *sig             = NULL;
  int  *signals_to_trap = NULL;
  int  *custom_sigs     = NULL;
  char *overrides_ptr   = NULL;
  char *overrides_copy  = NULL;
  char *token           = NULL;
  char *strtok_r_ptr    = NULL;

  const char delimeter[2] = ":";

  size_t num_sigs = 0;
  int i;

  if( inPar() == 1 )
  {
    /* Abort, because we are in a parallel region!
     * Some calls are not re-entrant
     */
    fprintf(stderr, "Error: signal_trap() called in a Parallel region!");
    signal_controlled_exit(1);
  }

  /* check whether the environment provides an override */
  overrides_ptr=getenv("UM_SIGNALS");

  if (overrides_ptr==NULL)
  {
    signals_to_trap=default_sigs;
  }
  else
  {
    /* override the default list */

    size_t overrides_len = 0;

    /* make a copy of the env var, so that we don't accidentally
     * change the environment when we call strtok_r
     */
    overrides_len = strlen(overrides_ptr);
    overrides_copy = malloc(overrides_len+1);
    strncpy(overrides_copy, overrides_ptr, overrides_len);

    if ( (sig_state.rank == 0) && (threadID() == 0) )
    {
      fprintf(stderr,
              "[%" PRId64 "] exceptions:"
              " Default handlers overriden with \"%s\"\n",
              sig_state.rank,overrides_copy);
    }

    /* Count how many tokens are in the string and
     *  allocate enough storage
     */
    token = strtok_r(overrides_copy, delimeter, &strtok_r_ptr);

    while( token != NULL )
    {
      num_sigs++;
      token = strtok_r(NULL, delimeter, &strtok_r_ptr);
    }

    custom_sigs = malloc(num_sigs+1);

    /* make a new copy, because strtok_r trashed the last one */
    strncpy(overrides_copy, overrides_ptr, overrides_len);

    /* add the tokens to the list to process */
    token = strtok_r(overrides_copy, delimeter, &strtok_r_ptr);
    num_sigs=0;

    while( token != NULL )
    {
      i = atoi(token);

      if( (i<1) || (i>31) )
      {
        fprintf(stderr,
                "[%" PRId64 "] exceptions:"
                " Signal \"%d\" is suspicious and ignored.\n",
                sig_state.rank, i);
      }
      else
      {
        if( (sig_state.rank==0) && (threadID()==0) )
        {
          fprintf(stderr,
                  "[%" PRId64 "] exceptions:"
                  " Override: Adding handler for signal %d.\n",
                  sig_state.rank, i);
        }

        custom_sigs[num_sigs] = i;
        num_sigs++;
      }

      token = strtok_r(NULL, delimeter, &strtok_r_ptr);
    }

    /* add a terminating zero */
    custom_sigs[num_sigs] = 0;

    if( (sig_state.rank==0) && (num_sigs==0)  && (threadID()==0) )
    {
      fprintf(stderr,
              "[%" PRId64 "] exceptions:"
              " Override: No signals are handled\n",
              sig_state.rank);
    }

    signals_to_trap=custom_sigs;
    free(overrides_copy);
  }

  sig_state.option=option;

  if( (sig_state.rank==0) && (sig_state.verbose) && (threadID()==0) )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " Setting option %" PRId64 "\n",
            sig_state.rank, option);
  }

  if( sig_state.option != signal_off )
  {
    memset(&sa, 0, sizeof(struct sigaction));

    sa.sa_flags = SA_SIGINFO;
    sa.sa_sigaction = signal_handler;
    sigemptyset(&sa.sa_mask);

    sig = signals_to_trap;

    while( *sig != 0 )
    {
      /* Add each signal to mask */
      sigaddset(&sa.sa_mask, *sig++);
    }

    sig = signals_to_trap;

    while( *sig != 0 )
    {
      if( *sig == SIGFPE )
      {
        extra_fp_options = true;
      }

      if( (sig_state.rank==0) && (sig_state.verbose) && (threadID()==0) )
      {
        fprintf(stderr,
                "[%" PRId64 "] exceptions:"
                " Registering signal %d with sigaction()\n",
                sig_state.rank, *sig);
      }

      sigaction(*sig++, &sa, NULL);
    }

    signal_trap_extra();

    gethostname(sig_state.hostname, MAX_HOSTNAME);
    sig_state.pid = getpid();
  }

  /* if we created a manual override, free the buffer */
  if( custom_sigs != NULL )
  {
    free(custom_sigs);
  }
}

/*----------------------------------------------------------------------------*/

/* Register a function to be called in a "bad situation" */
void signal_add_handler(void (*func)(void))
{
  signal_callback_t *cback;

  if (sig_state.verbose)
  {
    fprintf(stderr,"[%" PRId64 "] exceptions:"
            " Registering callback at 0x%08" PRIxPTR "\n",
            sig_state.rank,
            (uintptr_t)func);
  }

  if (sig_state.callbacks == NULL)
  {
    sig_state.callbacks = malloc(sizeof(*(sig_state.callbacks)));
    cback = sig_state.callbacks;
  }
  else
  {
    cback = sig_state.callbacks;

    while( cback->next != NULL )
    {
      cback = cback->next;
    }

    cback->next = malloc(sizeof(*(cback->next)));
    cback = cback->next;
  }

  cback->next = NULL;
  cback->callback = func;
}

/*----------------------------------------------------------------------------*/

/* Record a number to identify this task (usually an MPI rank) */
void signal_rank(int64_t pe)
{
  sig_state.rank = pe;
  threadFlush();
}

/*----------------------------------------------------------------------------*/

/* Enable debug output */
void signal_set_verbose(void)
{
  sig_state.verbose = 1;
}

/*----------------------------------------------------------------------------*/

/* Set the program name */
void signal_command(char *cmnd, int64_t len)
{
  if( len > MAX_FILENAME )
  {
    fprintf(stderr,
            "[%" PRId64 "] exceptions:"
            " Program name too long to store\n",
            sig_state.rank);
    return;
  }

  memcpy(sig_state.command, cmnd, (size_t)len);
  sig_state.command[len]='\0';
}

/*----------------------------------------------------------------------------*/

/* extra signal trap functions */
void signal_trap_extra(void)
{
  if (extra_fp_options)
  {
    /* set the floating point environment to the "FENV_ACCESS ON" state */
    FENVPRAGMA

    #if defined (IBM_XL_LIBC)

    signal_trap_extra_fp_ibm(&set_fpe_traps);

    #elif defined (GNU_LIBC)

    signal_trap_extra_fp_linux(&set_fpe_traps);

    #endif
  }
}
