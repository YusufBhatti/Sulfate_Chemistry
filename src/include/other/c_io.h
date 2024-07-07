#if !defined(C_IO_H) && defined(C95_2B)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* This Header is for Portio2B only */

#define C_IO_H
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/stat.h>
#include "pio_byteswap.h"
#include "portio_api.h"
#include "c_fort2c_types.h"

/*----------------------------------------------------------------------------*/
/* Enumerations and Typedefs                                                  */
/*----------------------------------------------------------------------------*/

/* enumerations */

enum roles         { c_io_user, c_io_readHelper , c_io_writeHelper, numRoles };

enum fileopenstate { oldFile, newFile, numOpenStates };

enum filemode      { readonly, readwrite, writeonly, numModes };

enum filestate     { idle, reading, writing, numStates };

enum outPutLevel   { outLevNone, outLevMin, outLevNormal, outLevOper,
                     outLevDiag, numOutLevels };

/* typedef enumerations */

typedef enum query_type {
  qry_file_exist
} query_type_t;

typedef enum HELPERSTATES {
  c_io_none,
  c_io_present,
  c_io_running,
  c_io_attaching,
  c_io_terminate,
  numHelperStates
} helperStates;

/* typedefs */

typedef struct c_io_layer_t {
  void  *parent;
  void  *next;
  void (* volatile helper_start)(void);
  void (* helper_stop)(void);
  volatile helperStates helper_state;
  pid_t   helper_tid;
  int64_t layer_number;
  char *name;
} c_io_layer_t;

typedef struct out_pack_bool_t {
  bool res;
} out_pack_bool_t;

typedef struct in_pack_fname_t {
  char *filename;
} in_pack_fname_t;

/* In order to avoid changing the interface, we will put all
   initialisation arguments in a struct.... */
typedef struct c_io_init_t {
  size_t       min_unit;
  size_t       max_unit;
  int64_t      wbuffering_buffer_size;
  int64_t      rbuffering_buffer_size;
  int64_t      rbuffering_buffer_count;
  bool         rbuffering_buffer_update;
  int64_t      rbuffering_buffer_prefetch;
  bool         timing_switch;
  c_io_layer_t *layers;
} c_io_init_t;

/* main attribute table */
typedef struct c_io_attributes_t {
  char *  filename;

  /* The file pointer for the unit */
  void *  fileHandle;

  /* either idle/reading/writing */
  int64_t mode;

  int64_t next_position;

  bool    modified_on_disk;

  /* Marker for IO events which modify the disk position to end-of-file.
   * Has a value of true if the last IO operation left the disk position marker
   * at the end of the file - e.g. an explicit seek to end; or a read which was
   * equal to or larger in size than the remaining data in the file.
   */
  bool    end_of_file;

} c_io_attributes_t;

/*----------------------------------------------------------------------------*/
/* Global variables and parameters                                            */
/*----------------------------------------------------------------------------*/

/* define common wait time for usleep() calls */

#define c_io_sleep_waittime 10

/* variables */

extern c_io_attributes_t *c_io_attributes;
extern c_io_init_t       *params;

extern int64_t c_io_attributes_size;
extern int64_t c_io_outLevel;

extern char modeText[numOpenStates][numModes][4];

extern bool c_io_initialised; /* is the library initialised? */

/* constants */

extern const int c_io_success;
extern const int c_io_fail;
extern const int c_io_file_end;
extern const int c_io_file_nopos;

/*----------------------------------------------------------------------------*/
/* Prototypes                                                                 */
/*----------------------------------------------------------------------------*/

/* DEPENDS ON: c_io.o */

/* Layer manipulations */

extern c_io_layer_t *c_io_get_layer        (int64_t);

extern c_io_layer_t *c_io_get_named_layer  (int64_t);

extern c_io_layer_t *c_io_get_bottom_layer (void);

extern c_io_layer_t *c_io_init_next_layer  (int64_t);

extern void          c_io_init_layer       (c_io_layer_t *,
                                            c_io_layer_t *,
                                            int64_t);

/* miscelany */

extern void          c_io_setOutPutLevel   (int64_t);

extern double        c_io_getTime          (void);

/* Layer zero API Prototypes */

/* File meta ops */

extern int64_t       c_io_open             (int64_t,
                                            char *,
                                            int64_t,
                                            int64_t);

extern int64_t       c_io_close            (int64_t);

extern int64_t       c_io_delete           (char *);

extern int64_t       c_io_sync             (int64_t);

extern int64_t       c_io_change_mode      (int64_t,
                                            int64_t);

extern int64_t       c_io_getExtent        (int64_t);

/* Actual IO */

extern int64_t       c_io_in               (int64_t,
                                            char *,
                                            int64_t,
                                            int64_t );

extern int64_t       c_io_out              (int64_t,
                                            char *,
                                            int64_t,
                                            int64_t );

extern int64_t       c_io_setpos           (int64_t,
                                            int64_t,
                                            int64_t);

extern int64_t       c_io_getpos           (int64_t,
                                            int64_t);

/*----------------------------------------------------------------------------*/
/* Convenience Functions                                                      */
/*----------------------------------------------------------------------------*/

/* IO bit-lengthed inline functions */

static inline int64_t c_io_in_8(int64_t unit, char *array, int64_t len)
{
  return c_io_in(unit, array, len, WORD8BYTES);
}

static inline int64_t c_io_in_32(int64_t unit, char *array, int64_t len)
{
  return c_io_in(unit, array, len, WORD32BYTES);
}

static inline int64_t c_io_in_64(int64_t unit, char *array, int64_t len)
{
  return c_io_in(unit, array, len, WORD64BYTES);
}

static inline int64_t c_io_out_8(int64_t unit, char *array, int64_t len)
{
  return c_io_out(unit, array, len, WORD8BYTES);
}

static inline int64_t c_io_out_32(int64_t unit, char *array, int64_t len)
{
  return c_io_out(unit, array, len, WORD32BYTES);
}

static inline int64_t c_io_out_64(int64_t unit, char *array, int64_t len)
{
  return c_io_out(unit, array, len, WORD64BYTES);
}

static inline int64_t c_io_setpos_8(int64_t unit, int64_t pos)
{
  return c_io_setpos(unit, pos, WORD8BYTES);
}

static inline int64_t c_io_setpos_32(int64_t unit, int64_t pos)
{
  return c_io_setpos(unit, pos, WORD32BYTES);
}

static inline int64_t c_io_setpos_64(int64_t unit, int64_t pos)
{
  return c_io_setpos(unit, pos, WORD64BYTES);
}

/* Other inline functions */

static inline int64_t c_io_isSuccessful(void)
{
  return c_io_success;
}

static inline int64_t c_io_isFailure(void)
{
  return c_io_fail;
}

static inline int64_t c_io_dummyInit(void * p)
{
  (void)(p);
  return c_io_success;
}

/* -------------------------------------------------------------------------- */

#endif

