#if !defined(PIO_BYTESWAP)
#define PIO_BYTESWAP

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* pio_byteswap.h */
/* DEPENDS ON: pio_byteswap.o */

#include <stdint.h>

/* -------------------------------------------------------------------------- */
/* Prototypes                                                                 */
/* -------------------------------------------------------------------------- */

typedef enum ENDIANNESS { bigEndian, littleEndian, numEndians } endianness;

extern int64_t pio_byteswap(void *array, int64_t len, int64_t word_len);
extern endianness get_machine_endianism(void);

#endif
