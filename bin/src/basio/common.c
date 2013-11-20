#if !defined(TYPES)
   #include "types.h"
#endif

#if !defined(COMMON)
   #define COMMON 1
#endif

// la fonction sleep():
#include <unistd.h> 

//-- globals ---------------------------------------------------------

int    Alloc_Borders(Bord_Idx n, Border **bb);
void   Flush_Segment(Segment *s);
int    Print_Segment(char *fname, Segment *s);
int    Input_Segment(char *fname, Segment *s);

void   Subtract_Counts(Count *n, Count *n1, Count *n2);
void   Add_Counts(Count *n, Count *n1, Count *n2);
void   Copy_Count(Count *source, Count *dst);
Count  Sum_Count(Count *c);
void   Fill_Counts(char *seq, Bord_Val k1, Count *c1, Bord_Val k2, Count *c2);

void   Copy_Border(Border *source, Border *dst);
int    Alphabet_Size(char *s);
int    Number_in_Alphabet(char c);
void   Change_Extension(char *in, char *out, char *ext);
int    Parse_Argv(int n, char *pp[], Comm_Line *c);
void   Print_Help(void);

Bord_Idx Count_Borders(FILE *fp);
Bord_Idx Distrib(Double t, Segment *s);

Int    Read_Sequence(char *fname, char **pp);

int 	 begin(char *s1, char *s2 );
long 	 Find_Seq_Start(FILE *fp);
Int    Count_Seq_Len(FILE *fp);

Int 	 Length_Between(Bord_Idx k1, Bord_Idx k2, Segment *s);


/*
Bord_Idx Min_Block_Len(Segment *s);
Bord_Idx Max_Block_Len(Segment *s);
Double   Lowest_Prob(Segment *s);
Double   Highest_Prob(Segment *s);
*/


void    Stat(Double *p, Int n, Statistics *s);
//int     Cmp_Double(Double *b1, Double *b2);
static int     Cmp_Double(void const *b1, void const *b2);
Double  Mean(Double *p, Int n);
Double  Stdev(Double *p, Int n);



//-- globals ---------------------------------------------------------

// definitions, not declarations;
// so other c-files can skip them if not needed

Alphabet Al;
int      Al_size;

//--------------------------------------------------------------------
int Alloc_Borders(Bord_Idx n, Border **bb)
{

   Bord_Idx i;

   if( !(*bb=(Border*)calloc(n, sizeof(Border))) ) {
      printf("\nError---Alloc_Borders(): alloc fail 1\n");
      //return NULL;
      return 0;
   } // end if

   for( i=0; i<n; i++) {
      if( !((*bb+i)->cnt=(Count*)calloc(Al_size, sizeof(Count))) ) {
         printf("\nError---Alloc_Borders(): alloc fail 2\n");
         //return NULL;
         return 0;
      } // end if
   } // end for

   return OK;

} // end Alloc_Borders()
//--------------------------------------------------------------------
void Flush_Segment(Segment *s)
{

   Bord_Idx i;

   for(i=0; i<s->num; i++) {
      free((s->b+i)->cnt);
   } // end for

   free(s->b);

   //s->num = NULL;
   s->num = 0;

} // end Flush_Segment()
//--------------------------------------------------------------------
// print to stdout if fname==NULL
int Print_Segment( char *fname, Segment *s)
{
   FILE *fp;
   Bord_Idx i;
   int      j;

//SC   out
//AL   a-c-g-t
//FN   0
//NR   6111
//VD
//BR        0    1.0000000000          0        0        0        0
//BR        3    0.7950589084          0        3        0        0
//BR        6    0.5906999406          0        4        1        1


   if( fname ) {
      if( !(fp=fopen(fname, "w")) ) {
         printf("\nError---Print_Segment(): can't open %s\n", fname);
         //return NULL;
         return 0;
      }
   } // end if(fname)
   else
      fp=stdout;

   fprintf( fp,

              "%s   %s"
            "\n%s   %s"
            "\n%s   %d"
            "\n%s   %d"
            "\n%s",

            SOURCE,   s->source,
            ALPHABET, s->a.p,
            FUNCTION, s->function,
            NUM_BORD, (int)(s->num),
            VOID

          );

   for( i=0; i<s->num; i++) {

      fprintf( fp, "\n%s %8d    %.*Lf  ", BORDER, (int)(s->b[i].k), DIGITS, s->b[i].prob );

      for( j=0; j < s->a.size; j++ ) {
         fprintf( fp, " %8d", (int)(s->b[i].cnt[j]) );
      } // end for j

   } // end for i

   if(fp!=stdout)
      fclose(fp);

   return OK;

} // end Print_Segment()
//--------------------------------------------------------------------
int Input_Segment( char *fname, Segment *s)
{
   FILE     *fp;
   int      j, k;

   Bord_Idx i=0, n;
   char     tmp[SSTR];

//SC   out
//AL   a-c-g-t
//FN   0
//NR   6111
//VD
//BR        0    1.0000000000          0        0        0        0
//BR        3    0.7950589084          0        3        0        0
//BR        6    0.5906999406          0        4        1        1


   if(!(fp=fopen(fname, "r"))) {
      printf("\nError---Input_Segment(): can't open %s\n", fname);
      //return NULL;
      return 0;
   } // end if

   if(fscanf(fp, "%*s%s%*s%s%*s%d", s->source, s->a.p, &(s->function) )!=3) {
      printf("\nError---Input_Segment(): invalid header in %s\n", fname);
      //return NULL;
      return 0;
   } // end if

   if((Al_size=(s->a.size)=Alphabet_Size(s->a.p))<2) {
      printf("\nError---Input_Segment(): invalid alphabet %s\n", s->a.p);
      //return NULL;
      return 0;
   } // end if

   if((n=Count_Borders(fp))<2) {
      printf("\nError---Input_Segment(): insufficient number of borders in %s\n", fname);
      //return NULL;
      return 0;
   } // end if

   rewind(fp);

   s->num = n;

   if( !Alloc_Borders(n, &(s->b)))
      //return NULL;
      return 0;

   while(fscanf(fp, "%s", tmp)==1) { // the one and only correct condition!!!

      if(strcmp(tmp, BORDER)) {
         continue;
      }
      else {

         k=0;

         k += fscanf(fp, "%d%f", (int*)&(s->b[i].k), (float*)&(s->b[i].prob) );

         for(j=0; j<Al_size; j++) {
            k += fscanf(fp, "%d", (int*)(s->b[i].cnt+j) );
         } // end for j

         if(k!=Al_size+2) {
            printf("\nError---Input_Segment(): invalid %s line\n", BORDER);
            //return NULL;
            return 0;
         } // end if

         i++;

      } // end else

   } // end while

   fclose(fp);

   strcpy(Al.p, s->a.p);

   return OK;

} // end Input_Segment()
//--------------------------------------------------------------------
// uses global: int Al_size
void Subtract_Counts(Count *n, Count *n1, Count *n2) // n=n1-n2
{
   int k;

   for(k=0; k<Al_size; k++)
      n[k]=n1[k]-n2[k];

} // end Subtract_Counts()
//--------------------------------------------------------------------
void Copy_Border(Border *source, Border *dst)
{

   dst->k     = source->k;
   dst->prob  = source->prob;

   Copy_Count(source->cnt, dst->cnt);

} // end Copy_Border()
//--------------------------------------------------------------------
void Copy_Count(Count *source, Count *dst)
{

   int i;

   for(i=0; i<Al_size; i++)
      dst[i] = source[i];

} // end Copy_Count()
//--------------------------------------------------------------------
// "a-c-g-t" -> 4, "ac-gt" -> 2, etc.
int Alphabet_Size(char *s)
{
   int j, n=1;

   for(j=0; s[j]; j++)
      if(s[j]==AL_SEPAR)
         n++;

   return n;

} // end Alphabet_Size()
//--------------------------------------------------------------------
int Parse_Argv(int n, char *pp[], Comm_Line *c ) // n==argc, pp==argv
{

   int  i,j;
   int  tot_non_opt=0;
   char *eq, *p, *opt;

   for(j=0; j<MAX_HASH; j++) // flush hash
      c->h[j].key=SPACE;

   // flush output
   c->infile[0]  = '\0';
   c->options[0] = '\0';

   opt = c->options;

   j=0;

   for(i=1; i<n; i++) {

      p = pp[i];

      if(p[0]=='-') { // this is an option

         if(eq=strchr(p,'=')) { // -t=float or --thresh=float

            c->h[j].key = tolower(p[(p[1]=='-')?2:1]);
            strcpy(c->h[j].value, eq+1);
            j++;

         } // end else if key-value
         else if(p[1]=='-') { // --long_form
            strncat(opt, p+2, 1);
         } // end --long_form
         else { // -p -a,  -pa --> "pa"
            strcat(opt, p+1);
         } // end else -pa

      } // end if option or -
      else { // not an option, apparently infile

         strcpy(c->infile, p);
         tot_non_opt++;

      } // end else

   } // end for

   //return (tot_non_opt==1)?OK:NULL; // error if not one non-option
   return (tot_non_opt==1)?OK:0; // error if not one non-option

} // end Parse_Argv()
//--------------------------------------------------------------------
// given counts for block (k0, k1), calculate counts for block (k0, k2),
// where k1<=k2. If c1==NULL, calculate counts for (k1, k2)
void Fill_Counts( char *seq, Bord_Val k1, Count *c1, Bord_Val k2, Count *c2 )
{

   Bord_Val j;
   int      n;

   if( c1 )
      Copy_Count(c1,c2);

   for( j=k1+1; j<=k2; j++) {
      if((n=Number_in_Alphabet(seq[j]))!=FAIL )
         c2[n]++;
   } // end for

}  // end Fill_Counts()
//--------------------------------------------------------------------
// uses global: Alphabet Al
int Number_in_Alphabet(char c)
{
   int i, k=0;

   for(i=0; i<strlen(Al.p); i++) {

      if(Al.p[i]==AL_SEPAR)
         k++;
      else if(tolower(c)==tolower(Al.p[i]))
         return k;

   } // end for

   printf("\nWarning---unknown letter %c", c);
   return FAIL;

} // end Number_in_Alphabet()
//--------------------------------------------------------------------
Int Read_Sequence(char *fname, char **pp)
{

   FILE   *fp;
   Int    seq_len, n;
   char   c, *p;
   long   seq_start;

   if( !(fp=fopen(fname, "r")) ) {
      printf("\nError---can't open sequence file %s\n", fname);
      //return NULL;
      return 0;
   } // end if

   if( (seq_start=Find_Seq_Start(fp))==-1 ) {
      printf("\nError---can't recognize file format\n");
      //return NULL;
      return 0;
   } // end if

   if(!(seq_len=Count_Seq_Len(fp))) {
      printf("\nError---sequence not found in %s\n", fname);
      //return NULL;
      return 0;
   } // end if


   if( !(p=(char*)calloc(seq_len+1, sizeof(char))) ) {
      printf("\nError---Read_Sequence(): alloc fail\n");
      //return NULL;
      return 0;
   } // end if

   n = 1; // sequence is in s[1], s[2], ..., s[n]

   fseek(fp, seq_start, SEEK_SET);

   while( (c = fgetc(fp))!= EOF ) {
      if( isalpha(c) ) {
         p[n++] = tolower(c);
      } // end if
   } // end while

   fclose(fp);

   *pp = p;

   return seq_len;

} // end Read_Sequence()
//--------------------------------------------------------------------
// copy in to out replacing (or adding) extension by ext
void Change_Extension(char *in, char *out, char *ext)
{

   char *dot;

   strcpy(out, in);

   if(dot=strrchr(out, '.'))
      *dot='\0';

   strcat(out, ".");
   strcat(out, ext);

} // end Change_Extension()
//--------------------------------------------------------------------
Count Sum_Count(Count *c)
{

   int   i;
   Count tot=0;

   for(i=0; i<Al_size; i++)
      tot += c[i];

   return tot;

} // end Sum_Count()
//--------------------------------------------------------------------
Bord_Idx Count_Borders(FILE *fp)
{

   Bord_Idx n=0;
   char     s[SSTR];

   while(fscanf(fp,"%s",s)==1) { // the only correct condition!!!
      if(!strcmp(s, BORDER))
         n++;
   } // end while

   return n;

} // end Count_Borders()
//--------------------------------------------------------------------
Bord_Idx Distrib(Double t, Segment *s)
{

   Bord_Idx i, n=0;

   for(i=0; i<s->num; i++) {
      if( s->b[i].prob<=t )
         n++;
   } // end for

   return n;

} // end Distrib()
//--------------------------------------------------------------------
// uses global: int Al_size
void Add_Counts(Count *n, Count *n1, Count *n2) // n=n1+n2
{
   int k;

   for(k=0; k<Al_size; k++)
      n[k]=n1[k]+n2[k];

} // end Add_Counts()
//--------------------------------------------------------------------
void Stat(Double *p, Int n, Statistics *s)
{

// Statistics: min, max, ave, var

   s->ave = Mean(p,n);
   s->dev = Stdev(p,n);

   //qsort((void*)p, n, sizeof(Double),(int(*)(void*,void*))Cmp_Double);
   qsort((void*)p, n, sizeof(Double), Cmp_Double);

   s->min = p[0];
   s->max = p[n-1];

} // end Statistics
//--------------------------------------------------------------------
//int Cmp_Double(Double *b1, Double *b2)
static int Cmp_Double(void const *p1, void const *p2)
{
	// +
	const Double *b1 = p1, *b2 = p2;
   if(*b1<*b2)
      return -1;
   else if(*b1>*b2)
      return 1;

   return 0;

} // end Cmp_Double()
//--------------------------------------------------------------------
Double Mean(Double *p, Int n)
{

   Double sum=0;
   Int    i;

   for(i=0; i<n; i++)
      sum += p[i];

   return sum/n;

} // end Mean()
//--------------------------------------------------------------------
// SUM(x-<x>)**2/(n-1) == (SUM(x**2)-n*<x>**2)/(n-1)
Double Stdev(Double *p, Int n)
{

   Double x, y, sum=0;
   Int    i;


   for(i=0; i<n; i++)
      sum += p[i]*p[i];

   x = Mean(p,n);

   y = (sum-n*x*x)/(n-1);

   if(y<=0)
      return 0;

   return sqrt(y);

} // end Stdev()
//--------------------------------------------------------------------
void Print_Help(void)
{

   printf(
#include  HELPFILE
         );

} // end Print_Help()
//--------------------------------------------------------------------
int begin(char *s1, char *s2 )
{
	return (strstr(s1, s2)==s1);
} // end begin();
//--------------------------------------------------------------------
long Find_Seq_Start(FILE *fp)
{

   char   s[LSTR];

   while( fgets(s, LSTR, fp) ) {
      if( begin(s, "SQ") || begin(s, "ORIGIN") || begin(s, ">") )
         return ftell(fp);
	} // end while

	return -1;

} // end Find_Seq_Start();
//--------------------------------------------------------------------
Int Count_Seq_Len(FILE *fp)
{

	//char c;
	int char_as_int;
	Int  n=0;

   //while((c=fgetc(fp))!= EOF) {
   while((char_as_int=fgetc(fp))!= EOF) {
	   char c = (char)char_as_int;
      if(isalpha(c)) {
         n++;
      } // end if
   } // end while

	return n;

} // end Count_Seq_Len()
//--------------------------------------------------------------------
Int Length_Between(Bord_Idx k1, Bord_Idx k2, Segment *s)
{
	return s->b[k2].k - s->b[k1].k;
} // end Length_Between()
//--------------------------------------------------------------------
