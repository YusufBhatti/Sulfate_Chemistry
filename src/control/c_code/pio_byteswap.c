/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* Description:                                                               */
/*   Portable byteswapping routines, which are used in Portio2A and Portio2B, */
/*   but which have minimal dependencies in order to aid use elsewhere. These */
/*   allow the checking of machine endianism, and accordingly, the changing   */
/*   of file endianism.                                                       */

#include <stddef.h>
#include "pio_byteswap.h"
#include "pio_byteswap_opt.h"

#if defined(_OPENMP) && defined(UM_USE_C_OPENMP_VIA_THREAD_UTILS)
#include "thread_utils.h"
#endif

/* We want the compiler to inline the par_swap??() family in the non-"OpenMP via
 * thread_utils" case for performance reasons.
 * However, in this case, the par_swap??() family must exist as seperate
 * routines in order to create the required function pointers.
 * Therefore, create an OpenMP and UM_USE_C_OPENMP_VIA_THREAD_UTILS dependent
 * macro to expand to "static inline" qualifiers only in the OpenMP via
 * thread_utils case.
 */

#if defined(_OPENMP) && defined(UM_USE_C_OPENMP_VIA_THREAD_UTILS)
#define INLINEQUAL
#else
#define INLINEQUAL static inline
#endif

/* -------------------------------------------------------------------------- */
/* Prototypes and Typedefs                                                    */
/* -------------------------------------------------------------------------- */

/* Prototypes */

/* Note: the first argument to these would more appropriately be (void *).
 * However, due to a bug, the Intel compiler does not recognise the VALUE
 * attribute in Fortran ABSTRACT INTERFACEs.
 * We must therefore use (void **) and pass-by-reference on the Fortran side.
 */
INLINEQUAL void par_swap64(void **const,
                           const int64_t *const restrict,
                           const int64_t *const restrict,
                           const int64_t *const restrict);

INLINEQUAL void par_swap32(void **const,
                           const int64_t *const restrict,
                           const int64_t *const restrict,
                           const int64_t *const restrict);

INLINEQUAL void par_swap16(void **const,
                           const int64_t *const restrict,
                           const int64_t *const restrict,
                           const int64_t *const restrict);

/* -------------------------------------------------------------------------- */
/* proceedure implentations                                                   */
/* -------------------------------------------------------------------------- */

int64_t pio_byteswap(void *array, int64_t len, int64_t word_len)
{
  int64_t status = 0;
  const int64_t imin = 0;
  const int64_t imax = len-1;
  const int64_t incr = 1;

#if defined(_OPENMP) && defined(UM_USE_C_OPENMP_VIA_THREAD_UTILS)
  /* intermediate restrict pointers are not needed for the
   * thread utils case.
   */
#else
  const int64_t *const restrict incr_ptr=&incr;
  const int64_t *const restrict imin_ptr=&imin;
  const int64_t *const restrict imax_ptr=&imax;
#endif

#if defined(_OPENMP) && !defined(UM_USE_C_OPENMP_VIA_THREAD_UTILS)
#if defined(_CRAYC)
  #pragma _CRI inline_always par_swap64, par_swap32, par_swap16
#endif
#endif

  if (word_len == 8)
  {
#if defined(_OPENMP) && defined(UM_USE_C_OPENMP_VIA_THREAD_UTILS)
    startOMPparallelfor(&array,par_swap64,&imin,&imax,&incr);
#else
    par_swap64(&array,imin_ptr,imax_ptr,incr_ptr);
#endif
  }
  else if (word_len == 4)
  {
#if defined(_OPENMP) && defined(UM_USE_C_OPENMP_VIA_THREAD_UTILS)
    startOMPparallelfor(&array,par_swap32,&imin,&imax,&incr);
#else
    par_swap32(&array,imin_ptr,imax_ptr,incr_ptr);
#endif
  }
  else if (word_len == 2)
  {
#if defined(_OPENMP) && defined(UM_USE_C_OPENMP_VIA_THREAD_UTILS)
    startOMPparallelfor(&array,par_swap16,&imin,&imax,&incr);
#else
    par_swap16(&array,imin_ptr,imax_ptr,incr_ptr);
#endif
  }
  else
  {
    status = 1;
  }

  return status;
}

/* -------------------------------------------------------------------------- */

INLINEQUAL void par_swap64(void **const array,
                           const int64_t * const restrict imin,
                           const int64_t * const restrict imax,
                           const int64_t * const restrict incr)
{
  {
    uint64_t *const restrict ptr_64 = (uint64_t *)*array;
    int64_t i;
    const int64_t span=(*imax-*imin)/(*incr);
    const int64_t min=*imin;
    const int64_t cincr=*incr;

#if defined(_OPENMP) && !defined(UM_USE_C_OPENMP_VIA_THREAD_UTILS)
    #pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<=span; i=i+1)
    {
      ptr_64[min+i*cincr] = (uint64_t)bswap_64(ptr_64[min+i*cincr]);
    }
  }
}

/* -------------------------------------------------------------------------- */

INLINEQUAL void par_swap32(void **const array,
                           const int64_t *const restrict imin,
                           const int64_t *const restrict imax,
                           const int64_t *const restrict incr)
{
  {
    uint32_t *const restrict ptr_32 = (uint32_t *)*array;
    int64_t i;
    const int64_t span=(*imax-*imin)/(*incr);
    const int64_t min=*imin;
    const int64_t cincr=*incr;

#if defined(_OPENMP) && !defined(UM_USE_C_OPENMP_VIA_THREAD_UTILS)
    #pragma omp parallel for default(shared) private(i)
#endif
    for (i=0; i<=span; i=i+1)
    {
      ptr_32[min+i*cincr] = (uint32_t)bswap_32(ptr_32[min+i*cincr]);
    }
  }
}

/* -------------------------------------------------------------------------- */

INLINEQUAL void par_swap16(void **const array,
                           const int64_t *const restrict imin,
                           const int64_t *const restrict imax,
                           const int64_t *const restrict incr)
{
  {
    uint16_t *const restrict ptr_16 = (uint16_t *)*array;
    int64_t i;
    const int64_t span=(*imax-*imin)/(*incr);
    const int64_t min=*imin;
    const int64_t cincr=*incr;

#if defined(_OPENMP) && !defined(UM_USE_C_OPENMP_VIA_THREAD_UTILS)
    #pragma omp parallel for default(none) private(i)
#endif
    for (i=0; i<=span; i=i+1)
    {
      ptr_16[min+i*cincr] = (uint16_t)bswap_16(ptr_16[min+i*cincr]);
    }
  }
}

/* -------------------------------------------------------------------------- */

endianness get_machine_endianism(void)
{
  endianness machine_endian;

  union {
    uint32_t thirtytwo_bit_int;
    uint8_t  four_byte_array[4];
  } byteorder={UINT32_C(0x01020304)};

  if (byteorder.four_byte_array[0]==1)
  {
    machine_endian=bigEndian;
  }
  else
  {
    machine_endian=littleEndian;
  }

  return machine_endian;
}

/* -------------------------------------------------------------------------- */
