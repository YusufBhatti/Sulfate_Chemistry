#if defined(C95_2B)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

#include "c_io_layers.h"

#define CIO_THIS_LAYER CIO_TRACE

#include "c_io_trace.h"

#include "c_io_nextlayer.h"

#include "c_io_internal.h"

#include <inttypes.h>

/* Functions in this file:
 *
 * Layer interfaces and functionality:
 *   init          : initialise structures and pass through
 *   fini          : report finialisation and pass through
 *   open          : report file details and pass through
 *   close         : report close operation and pass through
 *   in            : report bytes read and pass through
 *   out           : report bytes written and pass through
 *   getpos        : report position in words and pass through
 *   setpos        : report position in bytes and pass through
 *   change_mode   : report new mode and pass through
 *   sync          : report sync and pass through
 *   query         : report sync and pass through
 *
 * Implentation routines:
 *   [none]
 *
 * To make a new implementation of an IO layer, copy this
 * file, replacing lustreapi with your identifier.
 *
 */

/* -------------------------------------------------------------------------- */
/* Procedures                                                                 */
/* -------------------------------------------------------------------------- */

int64_t c_io_trace_init(c_io_init_t *init_params)
{
  c_io_printf(-1,"init: trace:");
  c_io_init_next_layer(CIO_NEXT_LAYER);

  return C_IO_INIT(init_params);
}

/* -------------------------------------------------------------------------- */

int64_t c_io_trace_fini(void)
{
  c_io_printf(-1,"fini: trace:");

  return C_IO_FINI();
}

/* -------------------------------------------------------------------------- */

int64_t c_io_trace_open(int64_t unit, char *filename, int64_t fileStatus,
                        int64_t fileMode)
{
  c_io_printf(unit,"open: trace: %s",filename);

  return C_IO_OPEN(unit,filename,fileStatus,fileMode);
}

/* -------------------------------------------------------------------------- */

int64_t c_io_trace_close(int64_t unit)
{
  c_io_printf(unit,"close: trace:");

  return C_IO_CLOSE(unit);
}

/* -------------------------------------------------------------------------- */

int64_t c_io_trace_in(int64_t unit, char *array, int64_t len, int64_t word_len)
{
  c_io_printf(unit,"in: trace: %" PRId64,(len*word_len));

  return C_IO_IN(unit,array,len,word_len);
}

/* -------------------------------------------------------------------------- */

int64_t c_io_trace_out(int64_t unit, char *array, int64_t len, int64_t word_len)
{
  c_io_printf(unit,"out: trace: %" PRId64,(len*word_len));

  return C_IO_OUT(unit,array,len,word_len);
}

/* -------------------------------------------------------------------------- */

int64_t c_io_trace_getpos(int64_t unit)
{
  int64_t r;

  r = C_IO_GETPOS(unit);
  c_io_printf(unit, "getpos: trace: %" PRId64, r);

  return r;
}

/* -------------------------------------------------------------------------- */

int64_t c_io_trace_setpos(int64_t unit, int64_t pos, int64_t word_len)
{
  c_io_printf(unit, "getpos: trace: %" PRId64, (pos*word_len));

  return C_IO_SETPOS(unit,pos,word_len);
}

/* -------------------------------------------------------------------------- */

int64_t c_io_trace_change_mode(int64_t unit, int64_t newMode)
{
  c_io_printf(unit, "change_mode: trace: %" PRId64, newMode);

  return C_IO_CHANGE_MODE(unit,newMode);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_trace_sync(int64_t unit)
{
  c_io_printf(unit, "sync: trace");

  return C_IO_SYNC(unit);
}

/*----------------------------------------------------------------------------*/

int64_t c_io_trace_query(int64_t unit,
                         query_type_t query,
                         void *in_pack,
                         void *out_pack)
{
  c_io_printf(unit, "query: trace: %" PRIdMAX, (intmax_t)query);

  return C_IO_QUERY(unit, query, in_pack, out_pack);
}

/*----------------------------------------------------------------------------*/

#endif
