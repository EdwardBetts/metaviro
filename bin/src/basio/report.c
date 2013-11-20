
#if !defined (TYPES)
   #include "types.h"
#endif

#if !defined(COMMON)
   #include "common.c"
#endif

#if !defined(MATHEM)
   #include "mathem.c"
#endif

//-- options and flags: --------------------------------------------

   static int     Format=0;
   static int     Print_Table=1;
   static int     Print_Summary=1;
   static int     Print_Blocks=0;

//-- globals: --------------------------------------------------------

   extern int   Al_size;
   Double       *Score;
   Double       Tot_score;
   char         *Seq;
   Int          Len;

   // no tabulation
   Int          Kmax=0;
   Double       *T=NULL;

   Int     MaxTab;
	int     TabStep; // tabulated are 0, TabStep, 2*TabStep, ...

//-- functions: ------------------------------------------------------
   int       Process_Options(Comm_Line *c);

   void     Print_Segm1(FILE *fp, Segment *s);
   void     Print_Segm2(FILE *fp, Segment *s);
   void     Print_Segm3(FILE *fp, Segment *s);
   void     Print_Segm4(FILE *fp, Segment *s);
   Double   Calculate_Scores(Segment *s, Double **p);

   void   Table(FILE *fp, Segment *s);
   void   Summary(FILE *fp, Segment *s);

	void 		Get_Freqs(Bord_Idx k1, Bord_Idx k2, Segment *s, Double *f);
	Double   Entropy(Double *f);



//--------------------------------------------------------------------
int main (int argc, char *argv[])
{

   Segment   s;
   FILE      *fp_out;
   char      outfile[LSTR];
   Comm_Line *cl;

   Seq=NULL;


   if(!(cl=(Comm_Line*)malloc(sizeof(Comm_Line)))) {
      printf("\nError---alloc fail\n");
      return EXIT_FAILURE;
   } // end if

   if(argc==1 || !Parse_Argv(argc,argv,cl) || !Process_Options(cl)) {
      Print_Help();
      return EXIT_SUCCESS;
   } // end if


   if(!Input_Segment(cl->infile, &s))
      return EXIT_FAILURE;

   Change_Extension(cl->infile, outfile, EXT2);
   if( !(fp_out=fopen(outfile, "w")) ) {
      printf("\nError---can't open %s\n", outfile);
      return EXIT_FAILURE;
   } // end if


   if(
         (Print_Blocks || Format) &&
         !(Len=Read_Sequence(s.source, &Seq))
     )
     return EXIT_FAILURE;

   if( !(Tot_score = Calculate_Scores(&s, &Score)))
      return EXIT_FAILURE;

   if(Print_Summary)
      Summary(fp_out, &s);

   if(Print_Table)
      Table(fp_out, &s);

   if(Format==1)
      Print_Segm1(fp_out, &s);
   else if(Format==2)
      Print_Segm2(fp_out, &s);
   else if(Format==3)
      Print_Segm3(fp_out, &s);
   else if(Format==4)
      Print_Segm4(fp_out, &s);

   if(fp_out!=stdout)
      fclose(fp_out);

   if(Seq)
      free(Seq);

   free(Score);

// free(LL);

   free(cl);
	return EXIT_SUCCESS;
} // end main()
//--------------------------------------------------------------------
void Summary(FILE *fp, Segment *s)
{

   Bord_Idx m;
   Bord_Val len;

// Statistics  stat_len;
// Double      *p_prob=NULL;
// Double      *p_len =NULL;

   m   = s->num-1;
   len = s->b[m].k-s->b[0].k;

/*
   if(m>1) {

      if(!(p_len=(Double*)calloc(m, sizeof(Double))) ) {
         printf("\nError---Summary(): alloc fail\n");
         return;
      } // end if

      for(i=0; i<m; i++)
         p_len[i] = s->b[i+1].k - s->b[i].k;

      Stat(p_len,m,&stat_len);

      free(p_len);

   } // end if n>1
   else {
      stat_len.min = stat_len.max = stat_len.ave = len;
      stat_len.dev = 0;
   } // end else
*/

/*



   if(!(p=(Double*)calloc(m, sizeof(Double))) ) {
      printf("\nError---Summary(): alloc fail\n");
      return;
   } // end if

*/

   fprintf(fp,

                 "%s Source    %s"
               "\n%s Alphabet  %s"
               "\n%s Length    %d"
               "\n%s"
               "\n%s Function  %d"
               "\n%s Score               %.6Lf"
               "\n%s Normalized score    %.6Lf"
               "\n%s"
//             "\n%s Block lengths       %d %d %d %f %f"
               "\n%s",

               SUMMARY, s->source,
               SUMMARY, s->a.p,
               SUMMARY, (int)len,
               SUMMARY,
               SUMMARY, s->function,
               SUMMARY, Tot_score,
               SUMMARY, Tot_score+(len)*log(Al_size),
               SUMMARY,
//             SUMMARY, m, (Int)stat_len.min, (Int)stat_len.max,
//                      stat_len.ave, stat_len.dev,
               VOID

          );


} // end Summary()
//--------------------------------------------------------------------
void Table(FILE *fp, Segment *s)
{

   Bord_Idx  i, j;
   int       k;
   Bord_Val  block_len;
	Double    *freqs;

//BLOCK       27 - 119      92       0.975562340600  1.163151 0.217391 0.163043 0.130435 0.489130

	if( !(freqs=(Double*)calloc(Al_size, sizeof(Double))) ) {
      printf("\nError---alloc fail\n");
      return;
	} // end if


   fprintf(fp, "\n%s     from   length    prob         score  entropy     -----%s contents------", HEADER, s->a.p);
   fprintf(fp, "\n%s", VOID);

   for(i=1; i<s->num; i++) {

      block_len = s->b[i].k-s->b[i-1].k;

		Get_Freqs(i-1,i,s, freqs);

      fprintf( fp,
                   "\n%s %8d    %-8d %.6Lf %9.3Lf   %2.3Lf  ",

                   BLOCK,
                   (int)s->b[i-1].k,
                   (int)block_len,
                   s->b[i].prob,
                   Score[i]+(block_len)*log(Al_size),
						 Entropy(freqs)

             );


      for( k=0; k<Al_size; k++ ) {
         fprintf( fp, "  %6.4Lf", (Double)(s->b[i].cnt[k] - s->b[i-1].cnt[k])/block_len );
      } // end for j

      if(Print_Blocks) {
         fprintf(fp, "  ");
         for(j=s->b[i-1].k+1; j<=s->b[i].k; j++)
            fprintf(fp, "%c", Seq[j]);
      } // end if

   } // end for

} // end Table()
//--------------------------------------------------------------------
void Print_Segm1(FILE *fp, Segment *s)
{

   Int    i, first, last, tot_printed=0;
   Border *p;
   int    j;

   first = s->b[0].k+1;       // number of first character to print
   last  = s->b[s->num-1].k;  // number of last  character to print

   p = s->b+1;

   fprintf(fp, "\n%s\n%s\n%s\n%s %5d ", VOID,VOID,VOID, SEQ, (int)first);

   for(i=first; i<=last; i++) {

      fprintf(fp, "%c", Seq[i]);
      tot_printed++;

      if( tot_printed%WIDTH==0 )
         fprintf(fp, " %-d\n%s %5d ", (int)i, SEQ, (int)i+1);

      if(i==p->k && i!=last) {

         fprintf(fp, "%c", SEQ_SEPAR);
         tot_printed++;

         if( tot_printed%WIDTH==0 )
            fprintf(fp, " %-d\n%s %5d ", (int)i, SEQ, (int)i+1);

         p++;

      } // end if

   } // end for

   // print last number
   for( j=0; j<WIDTH-tot_printed%WIDTH; j++)
      fprintf(fp, " ");

   fprintf(fp, " %-d\n%s", (int)last, VOID);

} // end Print_Segm1()
//--------------------------------------------------------------------
void Print_Segm2(FILE *fp, Segment *s)
{

   Bord_Idx  i;
   Int       j;


   fprintf(fp, "\n%s\n%s\n%s", VOID, VOID, VOID);

   for(i=1; i<s->num; i++) {
      fprintf(fp, "\n%s %5d--", SEQ, (int)s->b[i-1].k+1);
      for(j=s->b[i-1].k+1; j<=s->b[i].k; j++)
         fprintf(fp, "%c", Seq[j]);
      fprintf(fp, "--%-d ", (int)s->b[i].k);

   } // end for i

   fprintf(fp, "\n%s", VOID);


} // end Print_Segm2()
//--------------------------------------------------------------------
void Print_Segm3(FILE *fp, Segment *s)
{

   Bord_Idx n;
   int      *p;
   Int      start, x,i;

   n = s->b[s->num-1].k - s->b[0].k;
   start = s->b[0].k;

   if(!(p=(int*)calloc(n, sizeof(int)))) {
      printf("\nError---Print_Segm3(): alloc fail\n");
      return;
   } // end if

   for(i=1; i<s->num; i++)
      p[s->b[i].k-start]=1;

   fprintf(fp, "\n%s\n%s", VOID, VOID);


   for(x=0; x<n; x+=WIDTH) {

      fprintf(fp, "\n%s %5d ", SEQ, (int)(x+start+1));

      for(i=x+1; i<=x+WIDTH && i<=n; i++)
         fprintf(fp, "%c", Seq[i]);

      fprintf(fp, " %-d\n%s %5d ", (i==n+1)?(int)(start+n):(int)(x+start+WIDTH), SEQ, (int)(x+start+1));

      for(i=x+1; i<=x+WIDTH && i<=n; i++)
         fprintf(fp, "%c", p[i]?'|':'.');

      fprintf(fp, " %-d", (i==n+1)?(int)(start+n):(int)(x+start+WIDTH));

   } // end for

   fprintf(fp, "\n%s", VOID);

   free(p);

} // end Print_Segm3()
//--------------------------------------------------------------------
void Print_Segm4(FILE *fp, Segment *s)
{

   Bord_Idx i;
   Int      start, j;

   for(i=1; i<s->num; i++) {


      start = s->b[i-1].k;

      fprintf(fp, "\n%s\n%s %5d ", VOID, SEQ, (int)(start+1));

      for(j=start+1; j<=s->b[i].k; j++) {
         fprintf(fp, "%c", Seq[j]);
         if((j-start)%10==0) fprintf(fp, " ");
         if((j-start)%60==0) fprintf(fp, "%-d\n%s %5d ", (int)j, SEQ, (int)(j+1)); // end of string

      } // end for j

   } // end for i


} // end Print_Segm4()
//--------------------------------------------------------------------
Double Calculate_Scores(Segment *s, Double **pp)
{

   Double    (*f)(Count *c);
   Bord_Idx  i, n;
   Count     *c_tmp;
   Double    *p, tot;

   n = s->num;

   if( !(p=(Double*)calloc(n, sizeof(Double))) ) {
      printf("\nError---Calculate_Scores(): alloc fail 1\n");
      //return NULL;
      return 0;
   } // end if

   if(s->function==1)
      f = &f1;
   else
      f = &f2;

   if( !(c_tmp=(Count*)calloc(Al_size, sizeof(Count))) ) {
      printf("\nError---Calculate_Scores(): alloc fail 2\n");
      //return NULL;
      return 0;
   } // end if

   tot = 0;

   for(i=1; i<n; i++) {
      Subtract_Counts(c_tmp, (s->b+i)->cnt, (s->b+i-1)->cnt );
      tot += (p[i]=(*f)(c_tmp));
   } // end for

   *pp = p;
   return tot;

}  // end Calculate_Scores()
//--------------------------------------------------------------------
int  Process_Options(Comm_Line *c)
{
   int  i;
   char key;
   char val[SSTR];


   if(strchr(c->options,'h') || strchr(c->options,'?'))
      //return NULL;
      return 0;


   if(strchr(c->options,'t'))
      Print_Summary=0;

   if(strchr(c->options,'s'))
      Print_Table=0;

   if(strchr(c->options,'b'))
      Print_Blocks=1;


   for(i=0; (key=c->h[i].key)!=SPACE; i++) {

      strcpy(val, c->h[i].value);

      if(key=='f')
         Format = atoi(val);
      else
         //return NULL;
         return 0;

   } // end for

   return OK;

} // end Process_Options()
//--------------------------------------------------------------------
void Get_Freqs(Bord_Idx k1, Bord_Idx k2, Segment *s, Double *f)
{

	int k;
	Bord_Val len;

	len = s->b[k2].k-s->b[k1].k;

   for(k=0; k<Al_size; k++) {
      f[k] = (Double)(s->b[k2].cnt[k] - s->b[k1].cnt[k])/(s->b[k2].k - s->b[k1].k);
   } // end for k

} // end Get_Freqs()
//--------------------------------------------------------------------
Double   Entropy(Double *f)
{

   Double   theta, x=0;
	int   j;

   for( j=0; j<Al_size; j++)  {
		theta = f[j];
      x -= (theta>0)?theta*log(theta):0;
	} // end for

   return x;


} // end Entropy()
//--------------------------------------------------------------------
