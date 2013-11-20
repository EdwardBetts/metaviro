#if !defined(TYPES)
   #include "types.h"
#endif

//-----------------------------------------------------------------

   IntFloat  d2if(Double x);  // x -> ([x],{x})
   MantExp   d2me(Double x);
   Double    me2d(MantExp x);

   MantExp Sum_MantExp (MantExp x, MantExp y);
   MantExp Mult_MantExp(MantExp x, MantExp y);
   MantExp Exponent(Double x);

   void Print_MantExp(MantExp x);


//-----------------------------------------------------------------
IntFloat d2if(Double x)
{

   IntFloat a;

   a.i = (int)floor(x);
   a.f = x - a.i;

   return a;

} // end d2if()
//-----------------------------------------------------------------
MantExp d2me(Double x)
{

   IntFloat a;
   MantExp  y;

   a = d2if(log10(x));

   y.r = pow(10, a.f);
   y.c = a.i;

   return y;

} // end d2me()
//-----------------------------------------------------------------
Double me2d(MantExp x)
{
   return (x.r * pow(10, x.c));
} // end me2d()
//-----------------------------------------------------------------
MantExp Sum_MantExp (MantExp x, MantExp y)
{

   MantExp  a1, a2, z;
   IntFloat t;

   if(x.r==0)
      return y;
   else if(y.r==0)
      return x;

   a1 = (x.c>=y.c)?x:y;
   a2 = (x.c>=y.c)?y:x;

   t = d2if(log10(a1.r + a2.r*pow(10, a2.c-a1.c)));

   z.r = pow(10,t.f);
   z.c = a1.c + t.i;

   return z;

} // end Sum_MantExp()
//-----------------------------------------------------------------
MantExp Mult_MantExp(MantExp x, MantExp y)
{

   IntFloat a;
   MantExp  z;

   a = d2if( log10(x.r)+log10(y.r) );

   z.r = pow(10,a.f);
   z.c = x.c + y.c + a.i;

   return z;

} // end Mult_MantExp()
//-----------------------------------------------------------------
MantExp Exponent(Double x)
{

#define LOG10_E 0.434294481903 // log(e)

   MantExp   z;
   IntFloat  a;

   a = d2if(x*LOG10_E);

   z.r = pow(10,a.f);
   z.c = a.i;

   return z;


} // end Exponent()
//-----------------------------------------------------------------
void Print_MantExp(MantExp x)
{

   printf(">>> (%f %d)\n", (double)(x.r), (int)(x.c));


} // end ()
//-----------------------------------------------------------------
