#if !defined (TYPES)
   #include "types.h"
#endif

#if !defined(COMMON)
   #include "common.c"
#endif

#define MATHEM 1


//-- functions -------------------------------------------------------

   Double f1(Count *c);
   Double f2(Count *c);
   Double Info(Segment *s);

   Double SL(Int k1, Int k2);
   Double Sum_Log(Int k1, Int k2);
   int    Tabulate_SL(void);

// Double ge(Count *c);


//-- globals ---------------------------------------------------------

extern int     Al_size;

extern Double  *T;
extern Int     MaxTab; // max arg of Sum_Log
extern int     TabStep;

//--------------------------------------------------------------------
// tabulate Sum_log up to Kmax:
/*

MAX_MEM		maximum amount of memory allowed for tabulation of SumLog
MaxIdx		=MaxMem/sizeof(double), max allowed tabulation array size, sizeof(double)==8
MaxTab		array size required for tabulation
TabStep		best when ==1, since tabulated are TabStep, 2*TabStep, 3*TabStep, ... 

*/
int Tabulate_SL(void)
{

	Int i, LastIdx;
	Int MaxIdx = MAX_MEM/sizeof(Double)-1; // max allowed idx of Double array


	for( TabStep=1; MaxIdx*TabStep < MaxTab; TabStep++ )
		;

	LastIdx = MaxTab/TabStep;

	if(DEBUG) {
		printf("\nMAX_MEM=%d MaxIdx+1=%d MaxTab=%d TabStep=%d\n", MAX_MEM, (int)(MaxIdx+1), (int)(MaxTab), TabStep);
	} // end if DEBUG

	if(TabStep>1) {
		printf("\nWarning: TabStep==%d >1\n", TabStep);
	} // end if


   if(!(T=(Double*)calloc(LastIdx+1, sizeof(Double)))) {
      printf("\nError---Tabulate_SL(): alloc failed\n");
      //return NULL;
      return 0;
   } // end if

   T[0] = 0;

	if(TabStep==1) {  // save time on not calling SL()
	   for(i=1; i<=LastIdx; i++) {
	      T[i]=T[i-1]+log(i);
		} // end for
	} // end if
	else {
	   for(i=1; i<=LastIdx; i++) {
	      T[i]=T[i-1]+SL((i-1)*TabStep+1, i*TabStep);
		} // end for
	} // end else

   return OK;

} // end Tabulate_Sum_Log()
//--------------------------------------------------------------------
Double f1( Count *c )
{

   Double s=0;
   int    j;

   for(j=0; j<Al_size; j++)
      s += Sum_Log( c[j]+1, 2*c[j] );

   return
      s - Sum_Log( Sum_Count(c)+Al_size, 2*Sum_Count(c)+Al_size-1 );

} // end f1()
//--------------------------------------------------------------------
Double f2( Count *c )
{

   Double s=0;
   int    j;

   for(j=1; j<Al_size; j++)
      s += Sum_Log(1, c[j]);

   return
      s + Sum_Log(1, Al_size-1) - Sum_Log(c[0]+1, Sum_Count(c)+Al_size-1);

} // end f2()
//--------------------------------------------------------------------
/*
Double ge(Count *c) // Grassberger estimator
{
   Double a, s=0;
   int    i;

   for(i=0; i<Al_size; i++) {
      a  = (Double)(c[i]+1)/(Sum_Count(c)+Al_size);
      s -= a*log(a);
   } // end for

   return s;

} // end ge()
*/
//--------------------------------------------------------------------
Double SL(Int k1, Int k2)
{
   Int     i;
   Double  s=0;

   for(i=k1; i<=k2; i++)
      s+=log(i);

   return s;

} // end SL()
//--------------------------------------------------------------------
// Kmax==tabulation limit
Double Sum_Log(Int k1, Int k2)
{
	Int k = k1-1;

	// optimize esp. for tabstep==1

	if(TabStep==1)
		return T[k2]-T[k];

	return
		T[k2/TabStep] + SL(k2-k2%TabStep+1, k2) -
		(T[k/TabStep] + SL(k-k%TabStep+1, k));

/*

   if(k1>k2) {
      return 0;
   } // end if
   else {
      return
         (k1<Kmax?(T[k2>Kmax?Kmax:k2]-(k1<2?0:T[k1-1])):0) +
         (k2>=Kmax?SL(k1<Kmax?Kmax+1:k1,k2):0);
   } // end else
*/


} // end Sum_Log()
//--------------------------------------------------------------------
// information content of the whole segment
Double Info( Segment *s )
{

   Bord_Idx m;
   Bord_Val tot;
   int      j;
   Double   theta, x=0;


   m = s->num-1;
   tot = Sum_Count((s->b+m)->cnt);

   for( j=0; j<Al_size; j++) {

      theta = (Double)((s->b+m)->cnt[j])/tot;
      x -= (theta!=0)?theta*log(theta):0;

   } // end for

   return x;

} // end Info()
//--------------------------------------------------------------------
