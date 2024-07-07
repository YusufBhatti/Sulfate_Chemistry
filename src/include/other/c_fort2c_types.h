#if !defined(C_FORT2C_TYPES_H)
#define C_FORT2C_TYPES_H

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/* Declares variable type real to have the same wordlength          */
/* as a Fortran REAL data type, integer to have the same            */
/* wordlength as Fortran INTEGER data type.                         */

#include <stdint.h>
#include <inttypes.h>

#if defined(C_INT)

/* Fortran INTEGER is equivalent to C int */
typedef int integer;
typedef unsigned int u_integer;
#define PRIdinteger "d"

#elif defined(C_LONG_INT)

/* Fortran INTEGER is equivalent to C long int */
typedef unsigned long int u_integer;
typedef long int integer;
#define PRIdinteger "ld"

#elif defined(C_LONG_LONG_INT)

/* Fortran INTEGER is equivalent to C long long int */
typedef unsigned long long int u_integer;
typedef long long int integer;
#define PRIdinteger "lld"

#else

/* DEFAULT: Fortran INTEGER is equivalent to C int */
typedef unsigned int u_integer;
typedef int integer;
#define PRIdinteger "d"

#endif

/* We cannot tell what type of data is being written since we might have
 * packed data which is unsafe to be represented by a real datatype.  Using an
 * integer for um_data_t will allow us to not have to deal with unrepresentable
 * bit patterns in floating point numbers.
 */
#if defined(FRL8)
/* Fortran REAL is equivalent to C double */
typedef double real;
typedef int64_t um_data_t;
#else
/* Fortran REAL is equivalent to C float */
typedef float real;
typedef int32_t um_data_t;
#endif

/* Number of bytes in a word of a specific bit size */
#define WORD8BYTES  1 /*  8 bit */
#define WORD16BYTES 2 /* 16 bit */
#define WORD32BYTES 4 /* 32 bit */
#define WORD64BYTES 8 /* 64 bit */

#endif
