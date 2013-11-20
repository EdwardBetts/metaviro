
#if !defined (TYPES)
   #include "types.h"
#endif

#if !defined(COMMON)
   #include "common.c"
#endif

#if !defined(MATHEM)
   #include "mathem.c"
#endif

#include "mantexp.c"

//-- globals: ------------------------------------------------------

   extern int  Al_size;
   Double (*f)(Count *c);

   Double  *T;
   Int     MaxTab;
	int     TabStep; // tabulated are 0, TabStep, 2*TabStep, ...


//-- options and flags: --------------------------------------------

   static int     Function=2;
   static double  Bip=0;              // Border Insertion Penalty
   static int     Print_Ruler=0;
   static int     Show_Time=0;
   static int     Long_Range=0;
   static int     Calc_Segm=1;
   static int     Calc_Prob=1;
	static char    Outfile[SSTR];
	static Int     MBL=0;				// Minimal Block Length

//-- functions: ----------------------------------------------------

   int       Process_Options(Comm_Line *c);

   void      Start_Ruler(char *s);
   //void    Clear_Ruler(void);
   void      Update_Ruler(void);

   int       Find_pbi(Segment *s, Bord_Val *b_tmp);
   Bord_Idx  Condense_Borders(Segment *s, Bord_Val *p_in, Bord_Val **pp_out);
   int       Build_Segmentation(Segment *s);
   int       Calculate_Z(Segment *s);

//------------------------------------------------------------------
int main( int argc, char *argv[])
{

   Segment    s;
   Comm_Line *cl;
   Bord_Idx   old_num;
	Int        seq_len;


   if(!(cl=(Comm_Line*)malloc(sizeof(Comm_Line)))) {
      printf("\nError---alloc failed\n");
      return EXIT_FAILURE;
   } // end if

   if(argc==1 || !Parse_Argv(argc,argv,cl) || !Process_Options(cl)) {
      Print_Help();
      return EXIT_SUCCESS;
   } // end if

   if(Print_Ruler)
      setbuf(stdout, NULL);

   if(!Input_Segment(cl->infile, &s))
      return EXIT_FAILURE;

	// check command line option if function is not yet set
	if( s.function==0 ) {
		s.function=(Function==1)?1:2;
	} // end function not yet set

	f=(s.function==1)?&f1:&f2;

   old_num = s.num;

	seq_len = Length_Between(0, s.num-1, &s);

   if( seq_len>LONG_RANGE_THRESH) Long_Range =1;



   MaxTab = s.function==1?(2*seq_len+Al_size-1):(seq_len+Al_size-1);

   if(!Tabulate_SL()) {
      printf("\nError---tabulation failed\n");
      return EXIT_FAILURE;
	} // end if


   if(Calc_Segm && !Build_Segmentation(&s))
      return EXIT_FAILURE;

   if(Calc_Prob && !Calculate_Z(&s))
      return EXIT_FAILURE;


   Print_Segment(strlen(Outfile)?Outfile:cl->infile, &s);

   printf("BIP=%.2f, Borders: %d -> %d\n",
				Bip,
				(int)old_num,
				(int)s.num
			 );

   Flush_Segment(&s);

   if(Show_Time)
      printf("\nTime elapsed: %.2lf min\n", (double)clock()/CLOCKS_PER_SEC/SEC_IN_MIN );

   free(T);
   free(cl);

	return EXIT_SUCCESS;
} // end main()
//--------------------------------------------------------------------
int Build_Segmentation(Segment *s)
{

   Bord_Val  *tmp1, *tmp2;
   Bord_Idx  i, n_old, n_new;
   Border    *b_new;


   n_old = s->num;

   if( !(tmp1=(Bord_Val*)calloc(n_old, sizeof(Bord_Val))) ) {
      printf("\nError---Build_Segmentation(): alloc failed 1\n");
      //return NULL;
      return 0;
   } // end if

   Find_pbi(s, tmp1);

   // now tmp1 contains references to the best borders:
   //
   // -1 0 0 1 2 2 1 3 7 8   tmp1[i]
   //  0 1 2 3 4 5 6 7 8 9        i
   //
   // so the condensed tmp2 should be 0 1 3 7 8 9
   // (0 and 9 remain)

   n_new = Condense_Borders(s, tmp1, &tmp2);

   if(n_new < n_old) {  // replace border array

      if(!Alloc_Borders(n_new, &b_new))
         //return NULL;
         return 0;

      for( i=0; i<n_new; i++) {
         Copy_Border(s->b+tmp2[i], b_new+i);
      } // end for

      // copy tmp-segment to s
      Flush_Segment(s);
      s->num = n_new;
      s->b   = b_new;

   } // end if

   free(tmp1);
   free(tmp2);

   return OK;

} // end Build_Segmentation()
//-- globals: --------------------------------------------------------
unsigned long int   done, tot;
int  printed; // percents printed, incremented by 2
//--------------------------------------------------------------------
int Find_pbi(Segment *s, Bord_Val *b_tmp)
{

   Double    best_score, net_sc, sc_left, sc_right;
   Bord_Idx  best_pbi, n, i, j; // pbi == previous border index
   Count     *c_tmp;
   Double    *p;


   b_tmp[0] = 0;

   n = s->num;

   if(!(c_tmp=(Count*)calloc(Al_size, sizeof(Count))) ) {
      printf("\nError---Find_pbi(): alloc failed 1\n");
      //return NULL;
      return 0;
   } // end if

   if(!(p=(Double*)calloc(n, sizeof(Double))) ) {
      printf("\nError---Find_pbi(): alloc failed 2\n");
      //return NULL;
      return 0;
   } // end if

   tot  = (n-2)*(n-1)/2;

   if(Print_Ruler)
      Start_Ruler(s->source);

   for(i=1; i<n; i++) {

		if( MBL && Length_Between(0,i,s)<MBL ) {
			continue;
		} // end if

      // block as a whole
      best_score = (*f)((s->b+i)->cnt);
      best_pbi   = 0;

      for( j=1; j<i; j++ ) {

			if(MBL && (Length_Between(0,j,s)<MBL || Length_Between(j,i,s)<MBL)) {
				continue;
			} // end if

         sc_left  = p[j];

         // c_tmp = i-th count - j-th count
         Subtract_Counts( c_tmp, (s->b+i)->cnt, (s->b+j)->cnt );

         sc_right = (*f)(c_tmp);

         net_sc = sc_left + sc_right - Bip;

         if(net_sc > best_score) {
            best_score = net_sc;
            best_pbi   = j;
         } // end if

         done++; // number of cycles done
//printf("\ni=%-8d j=%-8d done=%-8d  tot=%-8d", i,j,done, tot);
//if(done%1000==0) printf("\ntot=%-8d done=%-8d", tot,done);

         if(Print_Ruler )
            Update_Ruler();

      } // end for( j=...

      p[i]     = best_score;
      b_tmp[i] = best_pbi;

//if(  i%10==0 ) printf("\ni=%-8d tot_i=%-8d", i,n-1);

   } // end for( i=...

   free(c_tmp);
   free(p);

   return OK;

} // end Find_pbi()
//--------------------------------------------------------------------
Bord_Idx Condense_Borders(Segment *s, Bord_Val *p_in, Bord_Val **pp_out)
{

   Bord_Idx n_in,n_out;
   Bord_Idx i,j;

   //  0 0 0 1 2 2 1 3 7 8   p_in, n_in==10
   //  0 1 2 3 4 5 6 7 8 9

   //  0 1 3 7 8 9           *pp_out, n_out==6
   //  0 1 2 3 4 5

   n_in  = s->num;
   n_out = 1;

   i=n_in-1;
   while( i ) {
      n_out++; // count remaining borders
      i=p_in[i];
   } // end while

   if( !(*pp_out=(Bord_Val*)calloc(n_out,sizeof(Bord_Val))) ) {
       printf("\nError---Condense_Borders(): can't allocate memory\n");
       return FAIL;
   } // end if

   (*pp_out)[n_out-1] = n_in-1; // maintain last border

   j = n_out-2;

   i=n_in-1;
   while( i ) {
      (*pp_out)[j--] = p_in[i];
      i=p_in[i];
   } // end while

   return n_out;

} // end Condense_Borders()
//--------------------------------------------------------------------
int Calculate_Z(Segment *s)
{


   Bord_Idx  i, j, m;
   Double    T, Inf;
   Bord_Val  kl, kr;
   Count     *c_tmp;
   Part_Func *z;
   MantExp   norm;


	m=s->num-1;


	if(m==1) { // nothing to calculate
		(s->b+0)->prob=(s->b+1)->prob=1;
		goto END;
	} // end if


   if(
         !(c_tmp=(Count*)calloc(Al_size, sizeof(Count)))    ||
         !(z=(Part_Func*)calloc(s->num, sizeof(Part_Func)))

     ) {
      printf("\nError---Calculate_Z(): alloc failed\n");
      //return NULL;
      return 0;
   } // end if


   if(!Long_Range)
      Inf = Info(s);

   if(Long_Range) {
      z[0].lr.r = z[m].rl.r = 1;
      z[0].lr.c = z[m].rl.c = 0;
   } // end if
   else {
      z[0].lr.r = z[m].rl.r = 1;
   } // end else

   tot = m*(m+1)/2;

   if(Print_Ruler)
      Start_Ruler(s->source);

   for( i=1; i<=m; i++ ) {
      for( j=0; j<i; j++ ) {

         kr = (s->b+i)->k;
         kl = (s->b+j)->k;

         if(!Long_Range)
            T = (kr-kl)*Inf-(i-j)*LN_2;

         Subtract_Counts(c_tmp, (s->b+i)->cnt, (s->b+j)->cnt);

         // to MantExp

         if(Long_Range) {
            z[i].lr = Sum_MantExp(z[i].lr, Mult_MantExp(z[j].lr,Exponent((*f)(c_tmp))));
         } // end if
         else {
            z[i].lr.r += z[j].lr.r * exp((*f)(c_tmp)+T);
         } // end else


         kr = (s->b+(m-j))->k;
         kl = (s->b+(m-i))->k;

         if(!Long_Range)
            T = (kr-kl)*Inf-(i-j)*LN_2;

         Subtract_Counts(c_tmp, (s->b+(m-j))->cnt, (s->b+(m-i))->cnt);

         // to MantExp

         if(Long_Range) {
            z[m-i].rl = Sum_MantExp(z[m-i].rl, Mult_MantExp(z[m-j].rl,Exponent((*f)(c_tmp))));
         } // end if
         else {
            z[m-i].rl.r += z[m-j].rl.r * exp((*f)(c_tmp)+T);
         } // end else

         done++; // number of cycles done

         if(Print_Ruler)
            Update_Ruler();


      } // end for j
   } // end for i

   if(Long_Range) { // similar to 1/z[m].lr.r
      norm.r = 1/(z[m].lr.r);
      norm.c = - z[m].lr.c;
   } // end if

   for( i=0; i<=m; i++ ) {

      if(Long_Range) {
         (s->b+i)->prob = me2d( Mult_MantExp(Mult_MantExp(z[i].lr,norm),z[i].rl ));
      } // end if
      else {
         (s->b+i)->prob =  z[i].lr.r / z[m].lr.r * z[i].rl.r;
      } // end else

   } // end for

   free(c_tmp);

END:

   if(Print_Ruler)
//    Clear_Ruler();
      printf("\n");

   return OK;

} // end Calculate_Z()
//--------------------------------------------------------------------
void Start_Ruler(char *s)
{
   done = 0;
   printed = 0;

   printf("\n%s %% 0", s);

} // end Start_Ruler()
//--------------------------------------------------------------------
/*
void Clear_Ruler(void)
{
   int i;

   printf("\r");

   for( i=0; i<LINE_WIDTH; i++ )
      printf(" ");

} // end Clear_Ruler()
*/
//--------------------------------------------------------------------
void Update_Ruler(void)
{

#define INCR 2

   int i;

   for( i=0; i<((double)done/tot*100-printed)/INCR; i++ ) {
      printf("%c", DOT);
      printed += INCR;
      if( printed%10==0 )
         printf("\b%d", printed); // \b: replace DOT by number
   } // end for

} // end Update_Ruler()
//--------------------------------------------------------------------
int  Process_Options(Comm_Line *c)
{
   int  i;
   char key;
   char val[SSTR];

	*Outfile='\0';

   if(strchr(c->options,'h') || strchr(c->options,'?'))
      //return NULL;
      return 0;


   if(strchr(c->options,'p'))
      Calc_Segm=0;

   if(strchr(c->options,'l'))
      Long_Range=1;

   if(strchr(c->options,'r'))
      Print_Ruler=1;

   if(strchr(c->options,'s'))
      Calc_Prob=0;

   if(strchr(c->options,'t'))
      Show_Time=1;

   for(i=0; (key=c->h[i].key)!=SPACE; i++) {

      strcpy(val, c->h[i].value);

      if(key=='b')
         Bip = atof(val);
      else if(key=='f')
         Function = atoi(val);
      else if(key=='m')
         MBL = atol(val);
      else if(key=='o')
         strcpy(Outfile, val);
      else
         //return NULL;
		 return 0;

   } // end for

   return OK;

} // end Process_Options()
//--------------------------------------------------------------------
