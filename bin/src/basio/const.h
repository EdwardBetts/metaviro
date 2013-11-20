
//--- portability ------------------------------------------------

#define THIS_IS_UNIX 1

#if !THIS_IS_UNIX
  #undef  EOF
  #define EOF 255
#endif

#define HELPFILE "_help_"


//--- constants --------------------------------------------------

#define DIGITS 12  // number of decimal digits printed in probabilities

#define FAIL          -1
#define OK             1
#define COMMENT       '#'
#define SPACE         ' '
#define DOT           '.'


#define AL_SEPAR      '-'
#define DEFAULT_AL    "a-c-g-t"
#define SEQ_SEPAR     '|'

#define LSTR   120
#define SSTR   80

// maximum amount of memory allowed for tabulation of SumLog
// sizeof(double)==8, so MAX_MEM should be at least 8 times
// greater than the length of a biggest sequence
#define MAX_MEM 600000000


//#define MAX_BORD      10000
#define LN_2          0.693147
#define SEC_IN_MIN    60
#define LINE_WIDTH    80

#define DEBUG         0

#define WIDTH         60

#define SGM_SUFFIX_WIDTH 5

#define LONG_RANGE_THRESH 9000


//--- fields, keywords and extensions ----------------------------

#define SEQ_START1  "SQ"
#define SEQ_START2  "ORIGIN"
#define EXT1       "sgm"
#define EXT2       "ant"

#define BORDER     "BR"
#define SOURCE     "SC"
#define ALPHABET   "AL"
#define FUNCTION   "FN"
#define NUM_BORD   "NR"
#define VOID       "VD"
#define SUMMARY    "SM"
#define HEADER     "HD"
#define BLOCK      "BL"
#define SEQ        "SQ"

//----------------------------------------------------------------
