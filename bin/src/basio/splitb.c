// 28.03.04: added `--shift' and `--part_length' option
// 30.03.04: fixed segment border error

#if !defined(TYPES)
   #include "types.h"
#endif

#if !defined(COMMON)
   #include "common.c"
#endif


//-- global variables: ---------------------------------------------

   extern Alphabet Al;
   extern int      Al_size;

//-- options and flags: --------------------------------------------

   static Int      Segment_Length=0;   // default is one segment
   static Int      Overlap=0;

   // Part_Offset and Part_Length make `split' work with substring of the original sequence
   static Int      Part_Offset=0;	// offset of the substring, value added to each border value
   static Int      Part_Length=0;	// length of the substring

   static int      Borders_After_Each_Lett=0;
   static char     Alt_Namespace[SSTR]; // Alternate namespace

//- functions: -----------------------------------------------------

   int   Process_Options(Comm_Line *c);
   void  Make_Segm_Fname(int counter, char *infile, char *outfile);
// int   Is_Int(char *s);
   int   Remove_Borders(Segment *s);

//------------------------------------------------------------------
int main (int argc, char *argv[])
{
	
   Comm_Line *cl;

   Segment  s;
   Bord_Val x, y, z;
   Bord_Idx d;
   Int      seqlen_raw, seqlen, i;
   int      segm_num=0;
   char     *seq_raw, *seq, outfile[LSTR];

   if(!(cl=(Comm_Line*)malloc(sizeof(Comm_Line)))) {
      printf("\nError---alloc fail\n");
      return EXIT_FAILURE;
   } // end if

   strcpy(Al.p, DEFAULT_AL); // initial setup
   Alt_Namespace[0] = '\0'; 

   if(argc==1 || !Parse_Argv(argc,argv,cl) || !Process_Options(cl)) {
      Print_Help();
      return EXIT_SUCCESS;
   } // end if

   if( (Al.size=Alphabet_Size(Al.p))<2 ) {
      printf("\nError---invalid alphabet %s\n", Al.p);
      return EXIT_FAILURE;
   } // end if

   Al_size = Al.size;

	// sequence is an array of chars, not a string terminated by '\0'
	if(!(seqlen_raw=Read_Sequence(cl->infile, &seq_raw)))
		return EXIT_FAILURE;

	// seqlen = min(Part_Length, seqlen_raw-Part_Offset);
	seqlen = !Part_Length || Part_Length>seqlen_raw-Part_Offset ? seqlen_raw-Part_Offset : Part_Length;

	if(!(seq=(char*)calloc(seqlen+1, sizeof(char)))) {
		printf("\nError---alloc fail\n");
		return EXIT_FAILURE;
	} // end if
	
	for(i=0; i<seqlen; i++)
		seq[i+1] = seq_raw[Part_Offset+i+1]; // sequence is in s[1], s[2]..., not in s[0], ...
	
	free(seq_raw);



	if(Segment_Length==0)
		Segment_Length=seqlen;


   //-- cycle over pieces ---------------------------------------
   for(x=y=0; x<seqlen && y<seqlen; x+=(Segment_Length-Overlap)) {    // without y<seqlen condition big overlaps work incorrectly

      y=x+Segment_Length;

      if(y>seqlen)
         y=seqlen; // truncate last block

      s.num = y-x+1;

      // this block has y-x+1 borders
      if(!Alloc_Borders(y-x+1, &(s.b)))
         return EXIT_FAILURE;

      s.function = 0;
      strcpy(s.a.p, Al.p);
      s.a.size = Al.size;
      strcpy(s.source, cl->infile);

      for(z=x; z<=y; z++) {   // z is assigned all border values

         d=z-x;               // d==shift from the first border, sort of border counter
         (s.b+d)->k = z;

         if(!d)               // don't calculate counts for first border,
            continue;         // since first count vector==0 after calloc

         Fill_Counts(seq, z-1, (s.b+d-1)->cnt, z, (s.b+d)->cnt);

      } // end for z

      if(Segment_Length!=seqlen) // more than one segment
         segm_num++;

      Make_Segm_Fname(segm_num, strlen(Alt_Namespace)?Alt_Namespace:cl->infile, outfile);

      if(!Borders_After_Each_Lett)
         Remove_Borders(&s);

	if(Part_Offset) 
		for(i=0; i<s.num; i++)
			s.b[i].k += Part_Offset; // same as (s.b+i)->k

      Print_Segment(outfile, &s);
      printf("%s: %d..%-d\t%d borders\n", outfile, (int)s.b[0].k, (int)s.b[s.num-1].k, (int)s.num);

      Flush_Segment(&s);

   } // end for x

   free(seq);
   free(cl);

   printf("\n");
   return EXIT_SUCCESS;

} // end main()
//--------------------------------------------------------------------
// aaa.bbb.seq -> aaa  -> aaa_0002.sgm
void Make_Segm_Fname(int counter, char *infile, char *outfile)
{
   char *x, s[LSTR];

   strcpy(s, infile);

   if(x=strchr(s,'.'))
      *x='\0';

   if(counter)
      sprintf(outfile, "%s_%0*d.%s", s, SGM_SUFFIX_WIDTH, counter, EXT1);
   else
      sprintf(outfile, "%s.%s", s, EXT1);

} // end Make_Segm_Fname()
//--------------------------------------------------------------------
/*
int Is_Int(char *s)
{

   char c;

   while( c=*s++ ) {
      if(!isdigit(c))
         return 0;
   } // end while

   return 1;

} // end Is_Int()
*/
//--------------------------------------------------------------------
int Remove_Borders(Segment *s)
{

   Border   *b_new;
   Bord_Idx i, j=0, old, new;
   int      k, n;
   Border   *p;

   old = s->num;

   new = 2;

   s->b->prob = (s->b+old-1)->prob = 1;

   for(i=1; i<old-1; i++) {
      p = s->b+i;
      n = 0;
      for(k=0; k<Al_size; k++) {
         n += ((p+1)->cnt[k]-p->cnt[k])*(p->cnt[k]-(p-1)->cnt[k]);
      } // end for

      if(!n) { // i-th border separates different characters
         new++;
         (s->b+i)->prob=1;
      } // end if

   } // end for

   if(!Alloc_Borders(new, &b_new))
      //return NULL;
      return 0;

   for(i=0; i<old; i++) {
      if((s->b+i)->prob) {
         (s->b+i)->prob = 0;
         Copy_Border(s->b+i, b_new+j++);
      } // end if
   } // end for

   Flush_Segment(s);
   s->num = new;
   s->b   = b_new;

   return new;

} // end Remove_Borders()
//--------------------------------------------------------------------
int  Process_Options(Comm_Line *c)
{
   int  i;
   char key;
   char val[SSTR];

   if(strchr(c->options,'h') || strchr(c->options,'?'))
      //return NULL;
      return 0;

   if(strchr(c->options, 'e'))
      Borders_After_Each_Lett=1;

   for(i=0; (key=c->h[i].key)!=SPACE; i++) {

      strcpy(val, c->h[i].value);

      if(key=='a')
         strcpy(Al.p, val);
      else if(key=='l')
         Segment_Length = atol(val);
      else if(key=='o')
         Overlap = atol(val);
      else if(key=='s')
         Part_Offset = atol(val);
      else if(key=='p')
         Part_Length = atol(val);
      else if(key=='n')
         strcpy(Alt_Namespace, val);
      else
         //return NULL;
		 return 0;

   } // end for

   return OK;

} // end Process_Options()
//--------------------------------------------------------------------
