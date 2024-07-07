#if defined(C95_2A) || defined(C95_2B)

/**********************************COPYRIGHT***********************************/
/*            (C) Crown copyright Met Office. All rights reserved.            */
/*         For further details please refer to the file COPYRIGHT.txt         */
/*        which you should have received as part of this distribution.        */
/**********************************COPYRIGHT***********************************/

/* Code Owner: Please refer to the UM file CodeOwners.txt                     */
/* This file belongs in section: C Code                                       */

/**********************************************************************/
/* Description:                                                       */
/*                                                                    */
/* Portable Data Conversion Routines                                  */
/*                                                                    */
/* Information:                                                       */
/*                                                                    */
/* This provides the Portable Data Conversion Routines:               */
/*                                                                    */
/*    IEEE2IBM  - Convert IEEE format reals/integers to IBM format    */
/*                                                                    */
/*    IBM2IEEE  - Convert IBM format reals/integers to IEEE format    */
/*                                                                    */
/*    IEEE2IEG  - Convert between IEEE formats for reals and integers */
/*                                                                    */
/*                                                                    */
/*    IEG2IEEE  - Convert between IEEE formats for reals and integers */
/*                                                                    */
/* It is intended that these routines be called as INTEGER functions  */
/* from Fortran.  Each returns an integer (non-zero) error code if    */
/* the conversion fails, with an error code of zero indicating a      */
/* successful conversion. See UMDP S5 for details of the error codes  */
/* returned and their meaning.                                        */
/*                                                                    */
/* These provide Fortran interfaces to the routines:                  */
/*                                                                    */
/*                                                                    */
/*    read_number  - read IEEE or IBM format number and return sign,  */
/*                   exponent and mantissa.                           */
/*    write_number - read sign, exponent and mantissa, make conversion*/
/*                   and output number in required IEEE or IBM format.*/
/*                                                                    */
/* Also provided are routines                                         */
/*                                                                    */
/*    MOVEBITS  - Routine to shift bits from one variable to another. */
/*                Provides functionality of CRAY routine MOVBITS.     */
/*                                                                    */
/*    MOVEBYTES - Routine to shift bytes from one variable to another.*/
/*                Provides functionality of CRAY routine STRMOV.      */
/*                                                                    */
/* Both MOVEBITS and MOVEBYTES are directly callable from Fortran.    */
/*                                                                    */
/* Code Description:                                                  */
/*   Language: C                                                      */
/*   This code is written to UMDP3 v6 programming standards.          */
/*                                                                    */
/**********************************************************************/

/* Standard header files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include <inttypes.h>

/* Include files for Fortran-to-C Interface                           */
#include "c_fort2c_prototypes.h"
#include "c_fort2c_types.h"
#include "pio_umprint.h"

/* Include file for Portable data conversion constants                */
#include "c_data_conv.h"

/* Prototype functions */
void read_number (dataformat format, size_t size, size_t offset,
                  uint64_t in_number, uint64_t *out_sign, int64_t *out_expo,
                  uint64_t *out_mant, size_t *mant_bits_read);

void write_number (dataformat format, size_t size, size_t offset,
                   uint64_t *out_number, uint64_t in_sign, int64_t in_expo,
                   uint64_t in_mant, size_t mant_bits_read);

#if defined(C_LOW)
void movebytes (integer src[],integer *isb,   integer *num,   integer dest[],
                integer *idb);
#elif defined(C_LOW_U)
void movebytes_ (integer src[],integer *isb,   integer *num,   integer dest[],
                 integer *idb);
#else
void MOVEBYTES (integer src[],integer *isb,   integer *num,   integer dest[],
                integer *idb);
#endif

#if defined(C_LOW)
void movebits (integer src[], integer *isb, integer *num,
               integer dest[],  integer *idb);
#elif defined(C_LOW_U)
void movebits_ (integer src[], integer *isb, integer *num,
                integer dest[],  integer *idb);
#else
void MOVEBITS (integer src[], integer *isb, integer *num,
               integer dest[],  integer *idb);
#endif

/**********************************************************************/
/* The following function is to be called from Fortran.  It           */
/* provides a portable version of Cray routine CRI2IEG -              */
/* interface to read_number/write_number                              */
/**********************************************************************/

int64_t ieee2ieg (datatypes *data_type, int64_t *num, void *ieg_num_out,
                  int64_t *offset_out, void *cri_num_in, int64_t *stride,
                  int64_t *size_num_in, int64_t *size_num_out)
{
/* Local variables */
int64_t errorcode=0 ;   /* Return error code (success=0) */
int i;
dataformat type_num_in, type_num_out;
int64_t offset;
uint64_t *cri_num_in_ptr, *ieg_num_out_ptr;

/* Variables in calls to read/write_number */
uint64_t sign=0;
int64_t expo=0;
uint64_t mant=0;
size_t mant_bits_read=0;

cri_num_in_ptr = (uint64_t *)cri_num_in;
ieg_num_out_ptr = (uint64_t *)ieg_num_out;

/* Check for valid input */

if (*num <= 0)
{
  char message[MAX_OUTSTR];

  errorcode = -3 ;
  snprintf(message, MAX_OUTSTR,
     "IEEE2IEG: Error - Invalid num = %" PRId64 ". Return code = %" PRId64,
     *num, errorcode);
  umPrint(message,"pio_data_conv");
  return errorcode ;
}

if (*stride != 1)
{
  char message[MAX_OUTSTR];

  errorcode = -7 ;
  snprintf(message, MAX_OUTSTR,
     "IEEE2IEG: Error - Invalid stride = %" PRId64 ". Return code = %" PRId64,
     *stride, errorcode);
  umPrint(message,"pio_data_conv");
  return errorcode ;
}

/* Convert Cray offset to my offset */

 offset = 8*WORD64BYTES - *offset_out - *size_num_out;

/* Check that offset is valid (must lie between 0 and word_length */
/* and be a multiple of the size of the output number             */

if (offset < 0 || offset%(*size_num_out) != 0
               || offset > 8*WORD64BYTES
               || *offset_out < 0)
{
  char message[MAX_OUTSTR];

  errorcode = -4 ;
  snprintf(message, MAX_OUTSTR,
     "IEEE2IEG: Error - Invalid bitoff = %" PRId64 ". Return code = %"
     PRId64, *offset_out, errorcode);
  umPrint(message,"pio_data_conv");
  return errorcode ;
}

/* Cray fortran uses IEEE data types */
/* Decide which data types are used */

if (*data_type == INTEGER){

  /* INTEGER data types */

  switch (*size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = IEEE16INT_FORMAT;  break;
    case 32:  type_num_in = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;
      snprintf(message, MAX_OUTSTR,
               "IEEE2IEG: Error - Invalid natlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode ;
    }
  }
  switch (*size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = IEEE16INT_FORMAT;  break;
    case 32:  type_num_out = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IEEE2IEG: Error - Invalid forlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else if (*data_type == REAL){

  /* REAL data types */

  switch (*size_num_in){               /* real IEEE in */
    case 32:  type_num_in = IEEE32FLT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64FLT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
               "IEEE2IEG: Error - Invalid natlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){              /* real IEEE out */
    case 32:  type_num_out = IEEE32FLT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64FLT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IEEE2IEG: Error - Invalid forlen %" PRId64 ". Return code = %"
               PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else if (*data_type == LOGICAL){

  /* LOGICAL data types can be treated as INTEGERs  */

  switch (*size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = IEEE16INT_FORMAT;  break;
    case 32:  type_num_in = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
               "IEEE2IEG: Error - Invalid natlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = IEEE16INT_FORMAT;  break;
    case 32:  type_num_out = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IEEE2IEG: Error - Invalid forlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else
{
  /* ERROR - unsupported type */
  char message[MAX_OUTSTR];

  errorcode = -2;

  snprintf(message, MAX_OUTSTR,
           "IEEE2IEG: Error - unsupported data type = %"
           PRIdMAX ". Return code = %" PRId64,
           (intmax_t)*data_type, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode;
}

/* Loop over all values converting them as we go along */
for (i=0; i < *num; i++){
  read_number(type_num_in, (size_t)*size_num_in, 0, *cri_num_in_ptr,
                     &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (size_t)*size_num_out, (size_t)offset,
               ieg_num_out_ptr, sign, expo, mant, mant_bits_read);

  if ( offset == 0 ) {
    offset = 8*WORD64BYTES - *size_num_out ;
    ieg_num_out_ptr++;
  }
  else {
    offset = offset - *size_num_out ;
  }

  cri_num_in_ptr++ ;

}

  return errorcode ;
}

/*********************************************************************/
/* The following function is to be called from Fortran.  It          */
/* provides a portable version of Cray routine IEG2CRI -             */
/* interface to read_number/write_number                             */
/*********************************************************************/

int64_t ieg2ieee (datatypes *data_type, int64_t *num, void *ieg_num_in,
                  int64_t *offset_in, void *cri_num_out, int64_t *stride,
                  int64_t *size_num_out, int64_t *size_num_in)
{
/* Local variables */
int64_t errorcode=0 ;  /* Return error code (success=0) */
int i;
dataformat type_num_in, type_num_out;
int64_t offset;
uint64_t *cri_num_out_ptr, *ieg_num_in_ptr;

/* Variables in calls to read/write_number */
uint64_t sign=0;
int64_t expo=0;
uint64_t mant=0;
size_t mant_bits_read=0;

cri_num_out_ptr = (uint64_t *)cri_num_out;
ieg_num_in_ptr = (uint64_t *)ieg_num_in;

/* Check for valid input */

if (*num <= 0)
{
  char message[MAX_OUTSTR];

  errorcode = -3;

  snprintf(message, MAX_OUTSTR,
           "IEG2IEEE: Error - Invalid num = %" PRId64 ". Return code = %"
           PRId64, *num, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode ;
}

if (*stride != 1)
{
  char message[MAX_OUTSTR];

  errorcode = -7;

  snprintf(message, MAX_OUTSTR,
           "IEG2IEEE: Error - Invalid stride = %" PRId64 ". Return code = %"
           PRId64, *stride, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode ;
}

/* Convert Cray offset to my offset */

 offset = 8*WORD64BYTES - *offset_in - *size_num_in;

/* Check that offset is valid (must lie between 0 and word_length */
/* and be a multiple of the size of the output number             */

if (offset < 0 || offset%( *size_num_in ) != 0
               || offset > 8*WORD64BYTES
               || *offset_in < 0)
{
  char message[MAX_OUTSTR];

  errorcode = -4;

  snprintf(message, MAX_OUTSTR,
           "IEG2IEEE: Error - Invalid bitoff = %" PRId64 ". Return code = %"
           PRId64, *offset_in, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode;
}


/* Cray fortran uses IEEE data types */
/* Decide which data types are used */

if (*data_type == INTEGER){

  /* INTEGER data types */

  switch (*size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = IEEE16INT_FORMAT;  break;
    case 32:  type_num_in = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
               "IEG2IEEE: Error - Invalid natlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = IEEE16INT_FORMAT;  break;
    case 32:  type_num_out = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IEG2IEEE: Error - Invalid forlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else if (*data_type == REAL){

  /* REAL data types  */

  switch (*size_num_in){               /* real IEEE in */
    case 32:  type_num_in = IEEE32FLT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64FLT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
               "IEG2IEEE: Error - Invalid natlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){              /* real IEEE out */
    case 32:  type_num_out = IEEE32FLT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64FLT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IEG2IEEE: Error - Invalid forlen %" PRId64 ". Return code = %"
               PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else if (*data_type == LOGICAL){

  /* LOGICAL data types can be treated as INTEGERs  */

  switch (*size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = IEEE16INT_FORMAT;  break;
    case 32:  type_num_in = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
               "IEG2IEEE: Error - Invalid natlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = IEEE16INT_FORMAT;  break;
    case 32:  type_num_out = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IEG2IEEE: Error - Invalid forlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else
{
  /* ERROR - unsupported data type */
  char message[MAX_OUTSTR];

  errorcode = -2;

  snprintf(message, MAX_OUTSTR,
           "IEG2IEEE: Error - unsupported data type = %"
           PRIdMAX ". Return code = %" PRId64,
           (intmax_t)*data_type, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode ;
}

/* Loop over all values converting them as we go along */

for (i=0; i<*num; i++){
  read_number(type_num_in, (size_t)*size_num_in, (size_t)offset,
              *ieg_num_in_ptr, &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (size_t)*size_num_out, 0, cri_num_out_ptr, sign,
               expo, mant, mant_bits_read);

  if ( offset == 0 ) {
    offset = 8*WORD64BYTES - *size_num_in ;
    ieg_num_in_ptr++ ;
  }
  else {
    offset = offset - *size_num_in ;
  }

  cri_num_out_ptr++ ;

}

return errorcode ;
}

/*********************************************************************/
/* The following function is to be called from Fortran.  It          */
/* provides a portable version of Cray routine CRI2IBM -             */
/* interface to read_number/write_number                             */
/*********************************************************************/

int64_t ieee2ibm (datatypes *data_type, int64_t *num, void *ibm_num_out,
                  int64_t *offset_out, void *cri_num_in, int64_t *stride,
                  int64_t *size_num_in, int64_t *size_num_out)
{
/* local variables */
int64_t errorcode=0 ;  /* Return error code (success=0) */
int i;
dataformat type_num_in, type_num_out;
int64_t offset;
uint64_t *cri_num_in_ptr, *ibm_num_out_ptr;

/* Variables in calls to read/write_number */
uint64_t sign=0;
int64_t expo=0;
uint64_t mant=0;
size_t mant_bits_read=0;

cri_num_in_ptr = (uint64_t *)cri_num_in;
ibm_num_out_ptr = (uint64_t *)ibm_num_out;

/* Check for valid input */

if (*num <= 0)
{
  char message[MAX_OUTSTR];

  errorcode = -3;

  snprintf(message, MAX_OUTSTR,
           "IEEE2IBM: Error - Invalid num = %" PRId64 ". Return code = %"
           PRId64, *num, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode ;
}

if (*stride != 1)
{
  char message[MAX_OUTSTR];

  errorcode = -7;

  snprintf(message, MAX_OUTSTR,
           "IEEE2IBM: Error - Invalid stride = %" PRId64 ". Return code = %"
           PRId64, *stride, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode ;
}

/* Convert Cray offset to my offset */
 offset = 8*WORD64BYTES - *offset_out - *size_num_out;

/* Check that offset is valid (must lie between 0 and word_length */
/* and be a multiple of the size of the output number             */

if (offset < 0 || offset%(*size_num_out ) != 0
               || offset > 8*WORD64BYTES
               || *offset_out < 0)
{
  char message[MAX_OUTSTR];

  errorcode = -4;

  snprintf(message, MAX_OUTSTR,
           "IEEE2IBM: Error - Invalid bitoff = %" PRId64 ". Return code = %"
           PRId64, *offset_out, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode;
}

/* Cray fortran uses IEEE data types */
/* decide which data types are used */

if (*data_type == INTEGER){

  /* INTEGER data types */

  switch (*size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = IEEE16INT_FORMAT;  break;
    case 32:  type_num_in = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
              "IEEE2IBM: Error - Invalid natlen = %" PRId64 ". Return code = %"
              PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = IEEE16INT_FORMAT;  break;
    case 32:  type_num_out = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IEEE2IBM: Error - Invalid forlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else if (*data_type == REAL){

  /* REAL data types */

  switch (*size_num_in){               /* real IEEE in */
    case 32:  type_num_in = IEEE32FLT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64FLT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
               "IEEE2IBM: Error - Invalid natlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){              /* real IBM out */
    case 32:  type_num_out = IBM32FLT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IBM64FLT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IEEE2IBM: Error - Invalid forlen %" PRId64 ". Return code = %"
               PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode ;
    }
  }
}
else if (*data_type == LOGICAL){

  /* LOGICAL data types can be treated as INTEGERs  */

  switch (*size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = IEEE16INT_FORMAT;  break;
    case 32:  type_num_in = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
               "IEEE2IBM: Error - Invalid natlen = %" PRId64
               " . Return code = %" PRId64,
               *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = IEEE16INT_FORMAT;  break;
    case 32:  type_num_out = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IEEE2IBM: Error - Invalid forlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_out,  errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else
{
  /* ERROR - unsupported type */
  char message[MAX_OUTSTR];

  errorcode = -2;

  snprintf(message, MAX_OUTSTR,
           "IEEE2IBM: Error - unsupported data type = %"
           PRIdMAX ". Return code = %" PRId64,
           (intmax_t)*data_type, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode;
}

/* Loop over all values converting them as we go along */

for (i=0; i < *num; i++){

  read_number(type_num_in,  (size_t)*size_num_in, 0, *cri_num_in_ptr,
                     &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (size_t)*size_num_out, (size_t)offset,
               ibm_num_out_ptr, sign, expo, mant, mant_bits_read);

  if ( offset == 0 ) {
    offset = 8*WORD64BYTES - *size_num_out ;
    ibm_num_out_ptr++ ;
  }
  else {
    offset = offset - *size_num_out ;
  }

  cri_num_in_ptr++ ;

}

  return errorcode ;
}

/*********************************************************************/
/* The following function is to be called from Fortran.  It          */
/* provides a portable version of Cray routine IBM2CRI -             */
/* interface to read_number/write_number                             */
/*********************************************************************/

int64_t ibm2ieee (datatypes *data_type, int64_t *num, void *ibm_num_in,
                  int64_t *offset_in, void *cri_num_out, int64_t *stride,
                  int64_t *size_num_out, int64_t *size_num_in)
{
/* Local variables */
int64_t errorcode=0 ;   /* Return error code (success=0) */
int i;
dataformat type_num_in, type_num_out;
int64_t offset;
uint64_t *ibm_num_in_ptr, *cri_num_out_ptr;

/* Variables in calls to read/write_number */
uint64_t sign=0;
int64_t expo=0;
uint64_t mant=0;
size_t mant_bits_read=0;

ibm_num_in_ptr = (uint64_t *)ibm_num_in;
cri_num_out_ptr = (uint64_t *)cri_num_out;

/* Check for valid input */

if (*num <= 0)
{
  char message[MAX_OUTSTR];

  errorcode = -3;

  snprintf(message, MAX_OUTSTR,
           "IBM2IEEE: Error - Invalid num = %" PRId64 ". Return code = %"
           PRId64, *num, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode;
}

if (*stride != 1)
{
  char message[MAX_OUTSTR];

  errorcode = -7;

  snprintf(message, MAX_OUTSTR,
           "IBM2IEEE: Error - Invalid stride = %" PRId64 ". Return code = %"
           PRId64, *stride, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode ;
}

/* Convert Cray offset to my offset */
 offset = 8*WORD64BYTES - *offset_in - *size_num_in;

/* Check that offset is valid (must lie between 0 and word_length */
/* and be a multiple of the size of the output number             */

if (offset < 0 || offset%(*size_num_in) != 0
               || offset > 8*WORD64BYTES
               || *offset_in < 0)
{
  char message[MAX_OUTSTR];

  errorcode = -4;

  snprintf(message, MAX_OUTSTR,
           "IBM2IEEE: Error - Invalid bitoff = %" PRId64 ". Return code = %"
           PRId64, *offset_in, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode ;
}

/* Cray fortran uses IEEE data types */
/* Decide which data types are used */

if (*data_type == INTEGER){

  /* INTEGER data types */

  switch (*size_num_in){                  /* IEEE integer in */
    case 16:  type_num_in = IEEE16INT_FORMAT;  break;
    case 32:  type_num_in = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
               "IBM2IEEE: Error - Invalid natlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){                 /* IEEE integer out */
    case 16:  type_num_out = IEEE16INT_FORMAT;  break;
    case 32:  type_num_out = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IBM2IEEE: Error - Invalid forlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else if (*data_type == REAL){

  /* REAL data types */

  switch (*size_num_in){                  /* real IBM in */
    case 32:  type_num_in = IBM32FLT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IBM64FLT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
              "IBM2IEEE: Error - Invalid natlen = %" PRId64 ". Return code = %"
              PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){                 /* real IEEE out */
    case 32:  type_num_out = IEEE32FLT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64FLT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
               "IBM2IEEE: Error - Invalid forlen %" PRId64 ". Return code = %"
               PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else if (*data_type == LOGICAL){

  /* LOGICAL data types can be treated as INTEGERs  */

  switch (*size_num_in){               /* IEEE integer in */
    case 16:  type_num_in = IEEE16INT_FORMAT;  break;
    case 32:  type_num_in = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_in = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -5;

      snprintf(message, MAX_OUTSTR,
               "IBM2IEEE: Error - Invalid natlen = %" PRId64 ". Return code = %"
               PRId64, *size_num_in, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
  switch (*size_num_out){              /* IEEE integer out */
    case 16:  type_num_out = IEEE16INT_FORMAT;  break;
    case 32:  type_num_out = IEEE32INT_FORMAT;  break;
#if defined(FRL8)
    case 64:  type_num_out = IEEE64INT_FORMAT;  break;
#endif
    default:
    {
      char message[MAX_OUTSTR];

      errorcode = -6;

      snprintf(message, MAX_OUTSTR,
              "IBM2IEEE: Error - Invalid forlen = %" PRId64 ". Return code = %"
              PRId64, *size_num_out, errorcode);
      umPrint(message,"pio_data_conv");

      return errorcode;
    }
  }
}
else
{
  /* ERROR - unsupported type */
  char message[MAX_OUTSTR];

  errorcode = -2;

  snprintf(message, MAX_OUTSTR,
           "IBM2IEEE: Error - unsupported data type = %"
           PRIdMAX ". Return code = %" PRId64,
           (intmax_t)*data_type, errorcode);
  umPrint(message,"pio_data_conv");

  return errorcode;
}

/* Loop over all values converting them as we go along */

for (i=0; i<*num; i++){
  read_number(type_num_in, (size_t) *size_num_in, (size_t) offset,
              *ibm_num_in_ptr, &sign, &expo, &mant, &mant_bits_read);
  write_number(type_num_out, (size_t) *size_num_out, 0, cri_num_out_ptr,
                      sign, expo, mant, mant_bits_read);

  if ( offset == 0 ) {
    offset = 8*WORD64BYTES - *size_num_in ;
    ibm_num_in_ptr++ ;
  }
  else {
    offset = offset - *size_num_in ;
  }

  cri_num_out_ptr++ ;

}
 return errorcode ;
}

/*********************************************************************/
/* The following function is to be called from C.  It                */
/* reads in numbers of different types:                              */
/* real/integer, IBM or IEEE  with offsets/packed etc                */
/* Called from the portable data conversion routines                 */
/*********************************************************************/

void read_number (dataformat format, size_t size, size_t offset,
                  uint64_t in_number, uint64_t *out_sign, int64_t *out_expo,
                  uint64_t *out_mant, size_t *mant_bits_read)
{

int64_t expo_bias=0;
uint64_t sign_mask=0, expo_mask=0, mant_mask=0, in_number_mask;
size_t i, first_bit;
size_t expo_bit_size=0;


switch (format){
  case IBM32FLT_FORMAT:
    /* format 1 IBM 32 bits */
    sign_mask =      ibm_sign_32;
    expo_mask =      ibm_expo_32;
    expo_bit_size =      ibm_expo_bits_32;
    expo_bias =      ibm_expo_bias_32;
    mant_mask =      ibm_mant_32;
    *mant_bits_read = ibm_mant_bits_32;
  break;

  case IBM64FLT_FORMAT:
    /* format 2 IBM 64 bits */
    sign_mask =      ibm_sign_64;
    expo_mask =      ibm_expo_64;
    expo_bit_size =      ibm_expo_bits_64;
    expo_bias =      ibm_expo_bias_64;
    mant_mask =      ibm_mant_64;
    *mant_bits_read = ibm_mant_bits_64;
  break;

  case IEEE32FLT_FORMAT:
    /* format 3  IEEE 32 bits */
    sign_mask =      ieee_sign_32;
    expo_mask =      ieee_expo_32;
    expo_bit_size =      ieee_expo_bits_32;
    expo_bias =      ieee_expo_bias_32;
    mant_mask =      ieee_mant_32;
    *mant_bits_read = ieee_mant_bits_32;
  break;

  case IEEE64FLT_FORMAT:
    /* format 4  IEEE 64 bits */
    sign_mask =      ieee_sign_64;
    expo_mask =      ieee_expo_64;
    expo_bit_size =      ieee_expo_bits_64;
    expo_bias =      ieee_expo_bias_64;
    mant_mask =      ieee_mant_64;
    *mant_bits_read = ieee_mant_bits_64;
  break;

  case IEEE16INT_FORMAT:
    /* format 6 IEEE 16 bit integers */
    sign_mask =      ieee_int_sign_16;
    mant_mask =      ieee_int_mant_16;
    *mant_bits_read = ieee_int_mant_bits_16;
    expo_bit_size = 0;
    expo_bias = 0;
  break;

  case IEEE32INT_FORMAT:
    /* format 7 IEEE 32 bit integers */
    sign_mask =      ieee_int_sign_32;
    mant_mask =      ieee_int_mant_32;
    *mant_bits_read = ieee_int_mant_bits_32;
    expo_bit_size = 0;
    expo_bias = 0;
  break;

  case IEEE64INT_FORMAT:
    /* format 8 IEEE 64 bit integers */
    sign_mask =      ieee_int_sign_64;
    mant_mask =      ieee_int_mant_64;
    *mant_bits_read = ieee_int_mant_bits_64;
    expo_bit_size = 0;
    expo_bias = 0;
  break;

  case FLT_TYPES:
    /* fall through */
  case INT_TYPES:
  break;
}

if ( (offset != 0)
     || (size != (*mant_bits_read + expo_bit_size + 1))
     || ((8*WORD64BYTES - offset - size) != 0) ){
  /* For offsets and size differences */

  /* shift number to right and clear off any other data in string */
  in_number = in_number >> offset;
  in_number_mask = ( (ONE64 << size) - 1 );
  in_number = in_number & in_number_mask;

  /* sign_mask */
  sign_mask = ONE64 << (size - 1);
}


/* For all types of nos */

/* Extract 1 bit sign and remove trailing 0s*/
  *out_sign= (in_number & sign_mask) >> (*mant_bits_read + expo_bit_size);

/* mantissa and expo */

if ((in_number & mant_mask) == 0 && (in_number & expo_mask) == 0) {
  *out_mant = 0 ;
  *out_expo = 0 ;
}
else{
if( format == IBM32FLT_FORMAT || format == IBM64FLT_FORMAT ){

  /* IBM type numbers */

  /*  find first non-zero bit
  (up to first three leading bits could be 0 from IBM)*/
  first_bit = 0;
  for (i = 1; first_bit == 0; i++){
    if((((in_number & mant_mask) >> (*mant_bits_read - i)) & 1 ) == 1){
      first_bit = i;
    }
    if(i == (*mant_bits_read)) first_bit = (*mant_bits_read);
  }

  /* shift mantissa to be of the form 1.fraction
     and scale expo appropriately */
  *out_mant = ((in_number & mant_mask) << first_bit);

  /*  IBM type - scale base 16 exponent to decimal and then correct
  for normalised mantissa */

  *out_expo = ( (( (int64_t)((in_number & expo_mask) >> *mant_bits_read)
                 - expo_bias) * 4) - (int64_t)first_bit);

  /* added extra bit to mantissa forming 1.frac */
  *mant_bits_read = *mant_bits_read + 1;

}
else if(format == IEEE32FLT_FORMAT || format == IEEE64FLT_FORMAT ){

  /* IEEE floating point number */

  /* For IEEE number, add in 1 above fractional mantissa */
  *out_mant = (in_number & mant_mask) | (ONE64 << *mant_bits_read);

  /* just expo, removes trailing 0s, removes offset*/
  *out_expo = ( (int64_t)((in_number & expo_mask) >> *mant_bits_read )
                - expo_bias);

  /* mant has grown by one bit */
  *mant_bits_read = *mant_bits_read + 1;
}
else {
/* For integers */
  *out_mant = mant_mask & in_number;
  *out_expo = 0;
}
}
/* NB - DATA FROM READ TO WRITE
sign:  1 bit representing sign of mantissa 0=+, 1=-
expo:  without any offset
mantissa: explicit 1 before fraction, normalised
mant_bits_read:  no. of bits in mantissa
offset:POSITIVE, no. of bits to left integer number is offset in string
size:  bit size of number read / wrote */

}

/*********************************************************************/
/* The following function is to be called from C.  It                */
/* write out numbers of different types:                             */
/* real/integer, IBM or IEEE with offsets/packed etc                 */
/* Called from the portable data conversion routines                 */
/*********************************************************************/

void write_number (dataformat format, size_t size, size_t offset,
                   uint64_t *out_number, uint64_t in_sign, int64_t in_expo,
                   uint64_t in_mant, size_t mant_bits_read)
{

  /* Local Variables used for the (bitwise) construction of output format */
  /* NB: must be exactly 64 bit                                           */
uint64_t sign_mask=0, expo_mask=0, mant_mask=0;
uint64_t sign=0, mant=0, conv_number, conv_number_mask;
uint64_t first_bit, temp;
int64_t expo_bias=0;
int64_t expo=0;

/* Other local variable                                                   */

size_t mant_bits_write=0, expo_bits_write=0;
size_t i;
size_t expo_diff;

switch(format){
  case IBM32FLT_FORMAT:
    /* format 1 IBM 32 bit */
    sign= in_sign << 31;
    mant_bits_write= ibm_mant_bits_32;
    expo_bits_write= ibm_expo_bits_32;
    sign_mask= ibm_sign_32;
    expo_mask= ibm_expo_32;
    mant_mask= ibm_mant_32;
    expo_bias= ibm_expo_bias_32;
  break;

  case IBM64FLT_FORMAT:
    /* format 2 IBM 64 bit */
    sign= in_sign << 63;
    mant_bits_write= ibm_mant_bits_64;
    expo_bits_write= ibm_expo_bits_64;
    sign_mask= ibm_sign_64;
    expo_mask= ibm_expo_64;
    mant_mask= ibm_mant_64;
    expo_bias= ibm_expo_bias_64;
  break;

  case IEEE32FLT_FORMAT:
    /* format 3 IEEE 32 bit */
    sign= in_sign << 31;
    mant_bits_write= ieee_mant_bits_32;
    expo_bits_write= ieee_expo_bits_32;
    sign_mask=ieee_sign_32;
    expo_mask= ieee_expo_32;
    mant_mask= ieee_mant_32;
    expo_bias= ieee_expo_bias_32;
  break;

  case IEEE64FLT_FORMAT:
    /* format 4 IEEE 64 bit */
    sign= in_sign << 63;
    mant_bits_write= ieee_mant_bits_64;
    expo_bits_write= ieee_expo_bits_64;
    sign_mask= ieee_sign_64;
    expo_mask= ieee_expo_64;
    mant_mask= ieee_mant_64;
    expo_bias= ieee_expo_bias_64;
  break;

  case IEEE16INT_FORMAT:
    /* format 6 IEEE 16 bit integers */
    sign= in_sign << ieee_int_mant_bits_16;
    mant_bits_write= ieee_int_mant_bits_16;
    expo_bits_write= 0;
    expo_mask= 0;
    expo_bias= 0;
    sign_mask= ieee_int_sign_16;
    mant_mask= ieee_int_mant_16;
  break;

  case IEEE32INT_FORMAT:
    /* format 7 IEEE 32 bit integers */
    sign= in_sign << ieee_int_mant_bits_32;
    mant_bits_write= ieee_int_mant_bits_32;
    expo_bits_write= 0;
    expo_mask= 0;
    expo_bias= 0;
    sign_mask= ieee_int_sign_32;
    mant_mask= ieee_int_mant_32;
  break;

  case IEEE64INT_FORMAT:
    /* format 8 IEEE 64 bit integers */
    sign= in_sign << ieee_int_mant_bits_64;
    mant_bits_write= ieee_int_mant_bits_64;
    expo_bits_write= 0;
    expo_mask= 0;
    expo_bias= 0;
    sign_mask= ieee_int_sign_64;
    mant_mask= ieee_int_mant_64;
  break;

  case FLT_TYPES:
    /* fall through */
  case INT_TYPES:
  break;
}

if ( in_mant == 0 && in_expo == 0 ) {
  mant = 0 ;
  expo = 0 ;
}
else {
/* write for real types - formats 1-4 inclusive */

if(format == IBM32FLT_FORMAT || format == IBM64FLT_FORMAT ){
  /* IBM numbers - reposition mantissa and re-scale expo to base 16 */

  /*  Find base 16 expo closest to but just larger than in_expo */

  /* expo_diff will always be in the range 1 <= expo_diff <= 4 */

  bool overflow=0;

  expo = in_expo/4 + (in_expo >= 0);

  expo_diff = (4 - ((uint64_t)in_expo & 3));
  // equivalent to expo_diff = (in_expo>=0)*4-(in_expo%4) on
  // a two's compliment system.
  // C99 standard requires two's compliment representation for fixed width
  // integral types (ie. int64_t).

  if ((expo_diff==4)&&(expo < 0))
  {
       expo++;
  }

  mant = in_mant >>  expo_diff;

  /* mant has shrunk back again, lose extra bit */

  mant_bits_read--;

  if ( expo >= expo_bias )
  {
    /* expo is too large for format                      */
    /* Set number to largest possible IBM floating point */
    expo = expo_bias - 1 ;
    mant = mant_mask ;
    overflow=1;
  }
  else if ( -expo > expo_bias )
  {
    /* expo < -expo_bias; ie. expo is too negative to fit format */
    /* Set number to smallest possible IBM floating point */
    expo = -(expo_bias);
    mant = 0 ;
    overflow=1;
  }

  expo = expo + expo_bias;

  if (!overflow){
  if (mant_bits_write < mant_bits_read)
  {
    mant = (mant >> (mant_bits_read - mant_bits_write - 1) ) ;

    if ( (mant & 1)  == 1 ) {

      if ( (mant >> 1) == mant_mask ) {
        /* Need to increment exponent and shrink mantissa accordingly */
        mant = (mant >> 1) + 1 ;
        mant = (mant >> 4) ;
        expo++ ;
      }
      else {
        mant = (mant >> 1) + 1 ;
      }
    }
    else {
      mant = (mant >> 1) ;
    }
  }
  else if(mant_bits_write > mant_bits_read)
  {
    mant = ( mant << (mant_bits_write - mant_bits_read) );
  }
  }
}
else if (format == IEEE32FLT_FORMAT || format == IEEE64FLT_FORMAT ){

  /* For IEEE numbers knock off the 1 before the decimal -
  mantissa is only the fraction and take one from mant_bits_read, as
  it effectively shrinks one bit */

  bool overflow=0;

  mant = in_mant ^ (ONE64 << (mant_bits_read - 1) );

  mant_bits_read--;

  if ( in_expo >= expo_bias ) {
    /* Set number to largest possible IEEE floating point */
    expo = 2*expo_bias + 1 ;  /* treat as in_expo = expo_bias + 1 */
    mant = mant_mask + 1 ;
    overflow = 1 ;
  }
  else if ( in_expo < -expo_bias ) {
    /* Set number to smallest possible IEEE floating point */
    expo = 0 ;               /* treat as in_expo = -expo_bias */
    mant = 0 ;
    overflow = 1;
  }
  else {
  /* Truncate/grow mantissa as appropriate to write_number format */
    expo = in_expo + expo_bias;
  }

  if (!overflow){
  if (mant_bits_write < mant_bits_read) {
    mant = (mant >> (mant_bits_read - mant_bits_write - 1) ) ;

    if ( (mant & 1)  == 1 ) {

      if ( (mant >> 1) == mant_mask ) {
        /* Need to increment exponent and shrink mantissa accordingly */
        mant = (mant >> 1) + 1 ;
        expo++ ;
      }
      else {
        mant = (mant >> 1) + 1 ;
      }
    }
    else {
      mant = (mant >> 1) ;
    }
  }
  else if(mant_bits_write > mant_bits_read)
  {
    mant = (mant << (mant_bits_write - mant_bits_read) );
  }
  }
}

else if(format>INT_TYPES){

  /* Write for integer types - formats 6-8 inclusive */
  /* NB integers are written in two's complement form.
     The Most Sig Bit contains sign info (0=+).
     Neg numbers are converse of their pos opposite + 1  */

  /* Positive integers */
  first_bit=0;
  if (in_sign == 0){
    /* Find first data bit */
    for (i=1;first_bit==0;i++){
      if ((( in_mant >> (mant_bits_read - i) ) & 1) == 1)
      {
        first_bit = i;
      }
      if (i == mant_bits_read) first_bit = mant_bits_read;
    }
  }

  /* Negetive integers */
  else{
    /* Find first data bit */
    temp = 0;
    for (i=1;first_bit==0;i++){
      if ( (in_mant >> (mant_bits_read - i)) < ( (temp<<1) + 1) )
      {
        first_bit = i;
      }
      temp = in_mant >> (mant_bits_read - i);
      if (i == mant_bits_read) first_bit = mant_bits_read;
    }
  }

  /* Shrink/truncate as appropriate */
  mant = (in_mant & mant_mask);

  /* If no is negative add ones above the Most Significant Bit
     to fill out two's complement form */

  /* Only need to fill out bits before bit 0 and bit mant_bits_read */

  if( (in_sign != 0)  &&  (mant_bits_write > mant_bits_read) ){
    mant = mant | ( ( (ONE64 << (mant_bits_write - mant_bits_read))
                     - 1 ) << mant_bits_read ) ;
    /* replace sign bit on mant */
    /* mant = mant | ( one << (mant_bits_read - 1) ); */
  }

  /* Warn if siginificant bits are lost
  mant_bits_diff = mant_bits_write - mant_bits_read;
  if ( (mant_bits_diff + first_bit) > 1) WARN
  */

  /* Null unused expo for build up number */
  expo = 0;
}
}

/* Build up converted number */
sign= sign_mask & sign;
mant= mant_mask & mant;
expo= (int64_t)expo_mask & (expo << mant_bits_write );
conv_number= sign | (uint64_t)expo | mant;

/* Dealing with different sizes and offsets within the longer
   data type */
/* i.e.,  format = 32 bit int, but size = 8, offset = 8 */
/*  ----------------11111111--------  */
/*  32 bit data     8b size  8b offest
    ^ this kind of packing is used in the UM */

  if ( (offset != 0) || (size != (mant_bits_write + expo_bits_write + 1))
      || ((8*WORD64BYTES - offset - size) != 0) ) {
  /* create mask for converted number */
    conv_number_mask = (ONE64 << size) - 1;
  }
  else {
    conv_number_mask = FULLMASK64 ;
  }

  /* apply offset to mask and number */
  conv_number_mask = conv_number_mask << offset;
  conv_number = conv_number << offset;

  /* write in the correct place within the original string */
  *out_number = conv_number | ( *out_number & (~conv_number_mask) );

}

/*********************************************************************/
/* The following function is to be called from Fortran.  It          */
/* provides a portable version of Cray routine STRMOV -              */
/* mirrored parameters and functionality: shifts bytes from          */
/* one variable to another.                                          */
/*********************************************************************/
void
#if defined(C_LOW)
movebytes
#elif defined(C_LOW_U)
movebytes_
#else
MOVEBYTES
#endif
(integer src[],  /* src:  source variable                       */
 integer *isb,   /* isb:  starting byte of moved byte,          */
                 /*       bytes numbered from 1 from left       */
 integer *num,   /* num:  number of bytes moved                 */
 integer dest[], /* dest: destianation variable                 */
 integer *idb)   /* idb:  starting byte in destination variable */
/* MOVEBYTES assumes kind=8 ie 64bit data */
{
  /* Local variables */
  integer  bytes, bytes_mask, dest_temp, one;
  integer  nbytes_to_move, num_left ;
  integer  src_start, dest_start ;
  integer  word_length=sizeof(integer) ; /* Length of word in bytes */
  int      i ;

  num_left = *num ;
  src_start = *isb ;
  dest_start = *idb ;

  /* Check which element we're reading from and adjust accordingly */
  for (src_start = *isb ; src_start>word_length ;
                          src_start=src_start-word_length) {
    src++ ;
  }

  /* Check which element we're writing to and adjust accordingly */
  for (dest_start = *idb ; dest_start>word_length ;
                           dest_start=dest_start-word_length) {
    dest++ ;
  }

  for (i=0 ; num_left > 0 ; i++) {

    /* First, set the nbytes_to_move to the number of bytes in
       src[src_index] between byte=src_start and byte=word_length */

    nbytes_to_move = word_length - src_start + 1;

    /* Check whether nbytes_to_move is greater than the number
       left to be moved.  If it is, then set equal to num_left */

    if (nbytes_to_move > num_left) nbytes_to_move = num_left ;

    /* Check whether there is enough room in dest[dest_index]
       between byte=dest_start and byte=8*word_length to insert
       nbytes_to_move and, if not, then set equal to the number
       of bytes available */

    if (nbytes_to_move > (word_length - dest_start + 1) ) {
       nbytes_to_move = word_length - dest_start + 1 ;
    }

    /* if (nbytes_to_move > (word_length - src_start + 1) ) {
       nbytes_to_move = word_length - src_start + 1;
    } */

    /* Move nbytes_to_move bytes from src_start in *src to dest_start
       in *dest  */

    /* make mask of same length as bits */

    if (nbytes_to_move == word_length) {
    /* This is effectively a copy dest[dest_index] = src[src_index] */
      bytes_mask = -1 ;
    }
    else{
      bytes_mask = 0;
      one = 1 ;
      bytes_mask = (one << (word_length * nbytes_to_move) ) - 1;
    }

    /* grab bits from src */

    bytes =
         (*src >>
             word_length*(8 - (src_start-1) - nbytes_to_move ) )
          & bytes_mask;

    /* put bits in dest */
    dest_temp = *dest  &
                ( ~(bytes_mask << word_length*( 8 -
                   (dest_start-1) - nbytes_to_move ) ) );

    *dest =  dest_temp |
        (bytes << word_length * ( 8 -
                    (dest_start-1) - nbytes_to_move ) );

    /* Update number of bytes moved and number left to move */
    num_left = num_left - nbytes_to_move ;

    /* Update src_start ready for move of next chunk and check
       whether we've moved onto the next data element  */

    src_start = src_start + nbytes_to_move ;

    if (src_start > word_length) {
      src_start = src_start - word_length ;
      /* Increment to next data element of src */
      src++ ;
    }

    /* Update dest_start ready for move of next chunk and check
       whether we've moved onto the next data element  */

    dest_start = dest_start + nbytes_to_move ;

    if (dest_start > word_length) {
      dest_start = dest_start - word_length ;
      /* Increment to next data element of dest */
      dest++ ;
    }
  }
}

/*********************************************************************/
/* The following function is to be called from Fortran.  It          */
/* provides a portable version of Cray routine MOVBIT -              */
/* mirrored parameters and functionality: shifts bits from           */
/* one variable to another.                                          */
/*********************************************************************/
void
#if defined(C_LOW)
movebits
#elif defined(C_LOW_U)
movebits_
#else
MOVEBITS
#endif
(integer src[],  /* src:  source variable                      */
 integer *isb,   /* isb:  starting bit of moved bits,          */
                 /*       bits numbered from 1 from left       */
 integer *num,   /* num:  number of bits moved                 */
 integer dest[], /* dest: destination variable                 */
 integer *idb)   /* idb:  starting bit in destination variable */
{
  /* Local variables */

  integer  bits, bits_mask, dest_temp, one;
  integer  nbits_to_move, num_left ;
  integer  src_start, dest_start ;
  integer  word_length=8*sizeof(integer) ; /* Word length in bits */
  int      i ;

  num_left = *num ;
  src_start = *isb ;
  dest_start = *idb ;

  /* Check which element we're reading from and adjust accordingly */
  for (src_start = *isb; src_start>word_length ;
                         src_start=src_start-word_length) {
    src++ ;
  }

  /* Check which element we're writing to and adjust accordingly */
  for (dest_start = *idb ; dest_start>word_length ;
                           dest_start=dest_start-word_length) {
    dest++ ;
  }

  for (i=0 ; num_left > 0 ; i++) {

    /* First, set the nbits_to_move to the number of bits in
       src[src_index] between bit=src_start and bit=word_length */

    nbits_to_move = word_length - src_start + 1;

    /* Check whether nbits_to_move is greater than the number
       left to be moved.  If it is, then set equal to num_left */

    if (nbits_to_move > num_left) nbits_to_move = num_left ;

    /* Check whether there is enough room in dest[dest_index]
       between bit=dest_start and bit=word_length to insert
       nbits_to_move and, if not, then set equal to the number
       of bits available */

    if (nbits_to_move > (word_length - dest_start + 1) ) {
       nbits_to_move = word_length - dest_start + 1 ;
    }

    /* Move nbits_to_move bits from src_start in *src to dest_start
       in *dest  */

    /* make mask of same length as bits */

    if (nbits_to_move == word_length) {
    /*  This is effectively a copy: dest[dest_index]=src[src_index] */
      bits_mask = -1 ;
    }
    else{
      bits_mask = 0;
      one = 1 ;
      bits_mask = (one << nbits_to_move) - 1;
    }

    /* grab bits from src */
    bits = (*src >> ( word_length - src_start - nbits_to_move + 1) )
         & bits_mask;

    /* put bits in dest */
    dest_temp = *dest
         & ( ~(bits_mask <<
              ( word_length - dest_start - nbits_to_move +1)) );
    *dest =  dest_temp |
                ( bits <<
                ( word_length - dest_start - nbits_to_move +1) );

    /* Update number of bits moved and number left to move */

    num_left = num_left - nbits_to_move ;

    /* Update src_start ready for move of next chunk and check
       whether we've moved onto the next data element  */

    src_start = src_start + nbits_to_move ;

    if (src_start > word_length) {
      src_start = src_start - word_length ;
      /* Increment to next data element of src */
      src++ ;
    }

    /* Update dest_start ready for move of next chunk and check
       whether we've moved onto the next data element  */

    dest_start = dest_start + nbits_to_move ;

    if (dest_start > word_length) {
      dest_start = dest_start - word_length ;
      /* Increment to next data element of dest */
      dest++ ;
    }
  }
}
#endif
