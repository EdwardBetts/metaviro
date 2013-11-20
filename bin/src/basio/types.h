
#define TYPES 1


#include <stdio.h>
#include <stdlib.h>  // calloc()
#include <ctype.h>   // isalpha(), etc.
#include <string.h>
#include <math.h>
#include <time.h>


#include "const.h"

//--------------- TYPE DEFINITIONS --------------------------------

typedef unsigned long int  Int;      // for element numbers
typedef long double        Double;

typedef Int   Count;      // for counts
typedef Int   Bord_Val;   // for border values
typedef Int   Bord_Idx;   // for border indexes

//--------------- STRUCTURES --------------------------------------
typedef struct {

   int   size;
   char  p[LSTR];  // "a-c-g-t" or "ac-gt"

} Alphabet;


typedef struct {

   Bord_Val  k;      // border
   Count     *cnt;   // count vector
   Double    prob;   // border probabilty

} Border;


typedef struct {

   Alphabet a;
   char     source[LSTR];   // source of sequence
   int      function;    // 0..2

   Bord_Idx num;         // length of array b
   Border   *b;          // b[0], b[1], ..., b[num-1]

} Segment;

//-----------------------------------------------------------------

typedef struct {

   Double min;
   Double max;
   Double ave;
   Double dev;

} Statistics;

//-----------------------------------------------------------------


typedef struct {

   long int  i;  // integer
   Double    f;  // remainder

} IntFloat; // Integer-Float: x = i + f = [x]+{x};


typedef struct {

   Double    r;    // mantissa
   long int  c;    // exponent

} MantExp; // (mantissa, exponent) representation of a long double


typedef struct {

   MantExp lr;
   MantExp rl;

} Part_Func; // (mantissa, exponent) representation of a long double


//-----------------------------------------------------------------

typedef struct {

   char key;
   char value[SSTR];

} Hash;

#define MAX_HASH 8

typedef struct {

   char     options[SSTR];  // -a -b -c, -abc, /a /b /c ... --> "abc"

   Hash     h[MAX_HASH];  // "-a=ac-gt"-->{'a',"ac-gt"}; "-t=0.95"-->{'t',"0.95"}

   char     infile[LSTR];

} Comm_Line;

//-----------------------------------------------------------------
