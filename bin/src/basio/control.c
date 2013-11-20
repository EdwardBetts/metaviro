#if !defined (TYPES)
   #include "types.h"
#endif

#if !defined(COMMON)
   #include "common.c"
#endif

//-- options and flags: --------------------------------------------

   static int    Print_Length=1;
   static int    Print_Prob=1;
   static Double Threshold=0;
   static int    Digits=DIGITS;
   static int    Verbose=0;


//-- functions: ----------------------------------------------------

   int  Process_Options(Comm_Line *c);

//--------------------------------------------------------------------
int main (int argc, char *argv[])
{

   Segment     s;
   Comm_Line   *cl;
   Statistics  stat;
   Double      *p=NULL;
   Bord_Idx    i, n;

   setbuf(stdout, NULL);

   if(!(cl=(Comm_Line*)malloc(sizeof(Comm_Line)))) {
      printf("\nError---alloc fail 1\n");
      return EXIT_FAILURE;
   } // end if

   if( !Parse_Argv(argc,argv,cl) || !Process_Options(cl)) {
      Print_Help();
      return EXIT_SUCCESS;
   } // end if


   if(!Input_Segment(cl->infile, &s))
      return EXIT_FAILURE;

   if(Threshold!=0)
     printf("\nBorders with P<=%f: %d\n", (double)Threshold, (int)Distrib(Threshold, &s));

   if(Print_Length || Print_Prob) {
      if(!(p=(Double*)calloc(s.num-1, sizeof(Double))) ) {
         printf("\nError---alloc fail 2\n");
         return EXIT_FAILURE;
      } // end if
   } // end if


   if(Print_Length) {

      if((n=s.num-1)>1) {

         for(i=0; i<n; i++)
            p[i] = s.b[i+1].k - s.b[i].k;

         Stat(p,n,&stat);

      } // end if n>1
      else {
         stat.min = stat.max = stat.ave = s.b[s.num-1].k - s.b[0].k;
         stat.dev = 0;
      } // end else

      printf(
               "Lengths:\t%s%d\t%s%d\t%s%d\t%s%.2Lf\t%s%.2Lf\n",

               (Verbose)?"num=":"",  (int)n,
               (Verbose)?"min=":"", (int)stat.min,
               (Verbose)?"max=":"", (int)stat.max,
               (Verbose)?"ave=":"",  stat.ave,
               (Verbose)?"dev=":"",  stat.dev
            );

   } // end if Print_Length

   if(Print_Prob) {

      if((n=s.num-2)>2) {

         for(i=0; i<n; i++)
            p[i] = s.b[i+1].prob;

         Stat(p,n,&stat);

      } // end if(n>2)
      else {
         stat.min = stat.max = stat.ave = 1;
         stat.dev = 0;
      } // end else

      printf( "Probabilities:");

      printf(  "\t%s%d\t%s%.*Lf\t%s%.*Lf\t%s%.*Lf\t%s%.*Lf",

				(Verbose)?"num=":"", (int)n,
				(Verbose)?"min=":"", Digits, stat.min,
				(Verbose)?"max=":"", Digits, stat.max,
				(Verbose)?"ave=":"", Digits, stat.ave,
				(Verbose)?"dev=":"", Digits, stat.dev
            );

   } // end if Print_Prob

   printf("\n");

   Flush_Segment(&s);

   if(p)
      free(p);

   free(cl);

	return EXIT_SUCCESS;
} // end main()
//--------------------------------------------------------------------
int  Process_Options(Comm_Line *c)
{

   int  i;
   char key;
   char val[SSTR];


   if(strchr(c->options,'h') || strchr(c->options,'?'))
      //return NULL;
      return 0;


   if(strchr(c->options, 'l'))
      Print_Prob=0;

   if(strchr(c->options, 'p'))
      Print_Length=0;

   if(strchr(c->options, 'v'))
      Verbose=1;

   for(i=0; (key=c->h[i].key)!=SPACE; i++) {

      strcpy(val, c->h[i].value);

      if(key=='t') {
         Threshold = atof(val);
         Print_Prob=Print_Length=0;
      }
      else if(key=='d')
         Digits = atoi(val);
      else
         //return NULL;
		 return 0;

   } // end for


   if(c->h[0].key=='t')
      Threshold = atof(c->h[0].value);

   return (Print_Length || Print_Prob || Threshold);

} // end Process_Options()
//--------------------------------------------------------------------
