#if !defined(EXCEPTIONS_H)
#define EXCEPTIONS_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* DEPENDS ON: exceptions.o */

#include <inttypes.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------*/
/* Prepocessor Defines                                                        */
/*----------------------------------------------------------------------------*/

/* we want as much hard coded and minimal use of dynamic arrays
 * in code expected to be called when the application is fragile
 */

/* Time to SIG_ALARM on an exception*/
#define ALARM_CALL 120
#define MAX_HOSTNAME 50
#define MAX_FILENAME 256
#define MAX_TRACEBACK 100

/*----------------------------------------------------------------------------*/
/* Enums & Typedefs                                                           */
/*----------------------------------------------------------------------------*/

/* Symbolic names for readability */
enum signal_options {signal_traceback, signal_core, signal_off};

/* types: -
 * signal_callback_t - a linked list of error handlers
 * signal_state_t    - a structure to hold all of our persistent state
 * fpetraps_t        - a structure of additional FP trap switches
 */

typedef struct SIGNAL_CALLBACK_T {
  void (*callback)(void);          /* a callback */
  struct SIGNAL_CALLBACK_T *next;  /* the next callback to do */
} signal_callback_t;


typedef struct SIGNAL_STATE_T {
  int64_t option;                  /* main option switch */
  int verbose;                     /* toggle noisy output */
  signal_callback_t *callbacks;    /* call back registrations */
  pid_t pid;                       /* the process id */
  int64_t rank;                    /* the MPI rank */
  char command[MAX_FILENAME];      /* the program */
  char hostname[MAX_HOSTNAME];     /* the host */
} signal_state_t;

typedef struct FPETRAPS_T {
  bool divbyzero;
  bool invalid;
  bool overflow;
  bool underflow;
  bool inexact;
} fpetraps_t;

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

/* variables */
extern signal_state_t sig_state;

/* functions */
extern int lineInfo(void *, char *, size_t, char *, size_t, int *);

/*----------------------------------------------------------------------------*/
/* Pragmas                                                                    */
/*----------------------------------------------------------------------------*/

/*
 * The define a macro for the C99 Standard pragma to modify the floating point
 * environment. (STDC FENV_ACCESS ON)
 */
#define FENVPRAGMA _Pragma("STDC FENV_ACCESS ON")

/* not all compilers support this pragma yet, so re-define it on problematic
 * ones as empty: -
 *   GCC, clang
 */

#if defined(__GNUC__) || defined(__clang__)
#undef FENVPRAGMA
#define FENVPRAGMA
#endif

#endif
