#if !defined (TYPES)
   #include "types.h"
#endif

#include "common.c"


//-- options and flags: --------------------------------------------

   // filtering criteria:

   static Double  Threshold=0;
   static Int     Number=0;
   static Double  Fraction=0;
	static char    Outfile[SSTR];

//- functions: -------------------------------------------------------

   int  Filter_Segment(Double t, Segment *s);
   int  Process_Options(Comm_Line *c);

//--------------------------------------------------------------------
int main (int argc, char *argv[])
{

   Segment   s;
   Comm_Line *cl;
   Bord_Idx  i, n, k;
   Double    *p;


   if(!(cl=(Comm_Line*)malloc(sizeof(Comm_Line)))) {
      printf("\nError---alloc fail\n");
      return EXIT_FAILURE;
   } // end if

   if(argc<3 || argc>4 || !Parse_Argv(argc,argv,cl) || !Process_Options(cl)) {
      Print_Help();
      return EXIT_SUCCESS;
   } // end if

   if(!Input_Segment(cl->infile, &s))
      return EXIT_FAILURE;

   n = s.num;

   if(Number || Fraction) { // convert Number or Fraction to Threshold

      k = (Number)?Number:floor(Fraction*n);

      if(k>n-2) {
//       printf("\nWarning---trying to filter out too much");
         k = n-2;
      } // end if n
      else if(k==0) {
//       printf("\nWarning---nothing filtered out");
         Flush_Segment(&s);
         return EXIT_FAILURE;
      } // end if n

      if(!(p=(Double*)calloc(n, sizeof(Double)))) {
         printf("\nError---alloc fail\n");
         return EXIT_FAILURE;
      } // end if

      for(i=0; i<n; i++) {
         p[i]=s.b[i].prob;
      } // end for

      qsort((void*)p, n, sizeof(Double), Cmp_Double);

      Threshold = p[k-1];

      free(p);

   } // end if(!Threshold)

   if(Filter_Segment(Threshold, &s))
   	Print_Segment(strlen(Outfile)?Outfile:cl->infile, &s);

   printf("\nBorders: %d -> %d\n", (int)n, (int)s.num);

   Flush_Segment(&s);

   free(cl);
	return EXIT_SUCCESS;
} // end main()
//--------------------------------------------------------------------
int Filter_Segment(Double t, Segment *s)
{

   Border   *b_new;
   Bord_Idx  i, j, old, new=0;

   old = s->num;

   // no probabilities assigned, nothing to filter out
   //if(s->b->prob==NULL)
   if(s->b->prob==0)
      return old;

   for(i=0; i<old; i++)
      if((s->b+i)->prob>t)
         new++;

   if(new<2 || new==old) {
//    printf("\nWarning---nothing filtered out for threshold %f\n", t);
      //return NULL;
      return 0;
   } // end if

   if(!Alloc_Borders(new, &b_new))
      //return NULL;
      return 0;

   j=0;

   for(i=0; i<old; i++) {
      if((s->b+i)->prob>t) {
         Copy_Border(s->b+i, b_new+j++);
      } // end if
   } // end for

   Flush_Segment(s);
   s->num = new;
   s->b   = b_new;

   return new;

} // end Filter_Segment()
//--------------------------------------------------------------------
int  Process_Options(Comm_Line *c)
{

   char key;
   char val[SSTR];
	int  i;

   if(strchr(c->options,'h') || strchr(c->options,'?'))
      //return NULL;
      return 0;

   for(i=0; (key=c->h[i].key)!=SPACE; i++) {

      strcpy(val, c->h[i].value);

	   if(key=='t')
	      Threshold = atof(val);
	   else if(key=='n')
	      Number = atol(val);
	   else if(key=='f')
	      Fraction = atof(val);
	   else if(key=='o')
	      strcpy(Outfile, val);
	   else
	      //return NULL;
		  return 0;

   } // end for

   return (Threshold || Number || Fraction);

} // end Process_Options()
//--------------------------------------------------------------------
