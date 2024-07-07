#if !defined(CIO_LAYERS_H)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */


#define CIO_LAYERS_H

/* These are the layers */

#include <stdint.h>
#include <sys/types.h>

#define CIO_NUMBERLAYERS 11

/* -------------------------------------------------------------------------- */
/*  Implemented "actual" layers                                               */
/* -------------------------------------------------------------------------- */

/* Actual layers have positive magic numbers. (The special CIO_NO_LAYER is 0).
 * Actual layers are defined in the order of the meta-layer they belong to
 * (higest meta-layer -> lowest meta-layer; see below for meta-layers).
 * When multiple actual layers are included for a meta-layer, they are included
 * in order of increasing value, as given here. For example, CIO_TIMING is
 * included before CIO_BYTESWAP (both have CIO_TOP_LAYER meta-layer) because the
 * value of CIO_TIMING (1) is less than the value of CIO_BYTESWAP (2). But both
 * are included before CIO_TRACE (3) because that belongs to the lower
 * CIO_INFO_LAYER meta-layer.
 */

#define CIO_NO_LAYER   0

#define CIO_TIMING     1

#define CIO_BYTESWAP   2

#define CIO_TRACE      3

#define CIO_WBUFFERING 4

#define CIO_RBUFFERING 5

#define CIO_BLACKHOLE  6

#define CIO_LUSTRE_API 7

#define CIO_THROTTLE   8

#define CIO_UNIX       9

#define CIO_LIBC       10

/* If we are the root layer, we neet to also define a table of strings
 * corresponding to the names of the layers
 */

#if defined(CIO_ROOT_LAYER)

static const char *c_io_layer_names[CIO_NUMBERLAYERS]={
  "NoLayer",
  "Timing",
  "ByteSwap",
  "Trace",
  "WBuffering",
  "RBuffering",
  "Blackhole",
  "LustreAPI",
  "Throttle",
  "POSIX",
  "LIBC"
};

#endif

/* -------------------------------------------------------------------------- */
/*  Meta-layers                                                               */
/* -------------------------------------------------------------------------- */

/* layers have positive magic numbers, so meta-layers must have negatives
 * These are defined in reverse order (lowest meta-layer -> higest meta-layer)
 * The special CIO_NO_METALAYER is always the lowest value.
 */

#define CIO_BOTTOM_LAYER -1
#define CIO_IOAPI_LAYER  -2
#define CIO_BUFFER_LAYER -3
#define CIO_INFO_LAYER   -4
#define CIO_TOP_LAYER    -5
#define CIO_NO_METALAYER -6

/* We also need a macro to return an actual layer for a given meta-layer's magic
 * number. The meta layers are assigned to actual layer placeholders. These are
 * set to one of the defined actual layers [see above] at the preprocessing
 * stage.
 *
 * The actual layer placeholders have the same name as the corresponding
 * meta-layer, but without the _LAYER suffix.
 */

#define META_LAYER_LOOKUP(metalayer) (                                         \
                                  metalayer == CIO_BOTTOM_LAYER ? CIO_BOTTOM : \
                                  metalayer == CIO_IOAPI_LAYER  ? CIO_IOAPI  : \
                                  metalayer == CIO_BUFFER_LAYER ? CIO_BUFFER : \
                                  metalayer == CIO_INFO_LAYER   ? CIO_INFO   : \
                                  metalayer == CIO_TOP_LAYER    ? CIO_TOP    : \
                                  CIO_NO_LAYER)

/* The layers have the following meta-layer assignments, which can be queried
 * with the GET_META_LAYER macro. Each layer should have a corresponding line
 * in the definition of GET_META_LAYER, of the form:
 *
 *     layer == <LAYER> ? <LAYER'S META-LAYER> : \
 *
 * This can be used later on to perform macro checks on the layer's
 * meta-layer, eg:
 *
 *     #if GET_META_LAYER(CIO_LIBC) == CIO_BOTTOM_LAYER
 *
 * which would evaluate to true in this example.
 */

#define GET_META_LAYER(layer) (   layer == CIO_TIMING     ? CIO_TOP_LAYER    : \
                                  layer == CIO_BYTESWAP   ? CIO_TOP_LAYER    : \
                                  layer == CIO_TRACE      ? CIO_INFO_LAYER   : \
                                  layer == CIO_WBUFFERING ? CIO_BUFFER_LAYER : \
                                  layer == CIO_RBUFFERING ? CIO_BUFFER_LAYER : \
                                  layer == CIO_UNIX       ? CIO_BOTTOM_LAYER : \
                                  layer == CIO_LIBC       ? CIO_BOTTOM_LAYER : \
                                  layer == CIO_BLACKHOLE  ? CIO_IOAPI_LAYER  : \
                                  layer == CIO_THROTTLE   ? CIO_IOAPI_LAYER  : \
                                  layer == CIO_LUSTRE_API ? CIO_IOAPI_LAYER  : \
                                  CIO_NO_METALAYER )

#endif
