#if !defined(TYPES)
   #include "types.h"
#endif

#if !defined(COMMON)
   #include "common.c"
#endif

static int Cmp_Bord_Val(void const *b1, void const *b2);
int Process_Options(Comm_Line *c);

int Fast_Mode=1;

//--------------------------------------------------------------------
int main (int argc, char *argv[])
{

	int counter=1;
	char *seq, file[SSTR], tmp[SSTR], outfile[SSTR];
	Comm_Line *cl;
	Bord_Idx  i, j, num_bord, tot_bord=0;
	Bord_Idx	db; // distinct borders;
	FILE *fp;
	Bord_Val *p, *p0;
	Segment s_tmp, s;

	if(!(cl=(Comm_Line*)malloc(sizeof(Comm_Line)))) {
		printf("\nError---alloc fail\n");
		return EXIT_FAILURE;
	} // end if

	if(argc==1 || !Parse_Argv(argc,argv,cl) || !Process_Options(cl)) {
		Print_Help();
		return EXIT_SUCCESS;
	} // end if

	// count the borders:
	while(1) {
		// cl->infile is a namespace here
		sprintf(file, "%s_%0*d.%s", cl->infile, SGM_SUFFIX_WIDTH, counter, EXT1);
	   	if(!(fp=fopen(file, "r")))
			break;
		if( Fast_Mode ) { // trust the NR-line
			fscanf(fp, "%*s%*s%*s%*s%*s%*s%*s%d", &num_bord);
		} // end if
		else {
			num_bord = Count_Borders(fp);
		} // end else
		tot_bord += num_bord;
		counter++;
		fclose(fp);
	} // end while


	// read the borders:

	if( !(p=(Bord_Val*)calloc(tot_bord+1, sizeof(Bord_Val))) ) {
		printf("\nError---main(): alloc fail\n");
		return EXIT_FAILURE;
	} // end if

	p0 = p;

	counter=1;

	while(1) {

		// cl->infile is a namespace here
		sprintf(file, "%s_%0*d.%s", cl->infile, SGM_SUFFIX_WIDTH, counter, EXT1);

	   	if(!(fp=fopen(file, "r")))
				break;

	   	while(fscanf(fp, "%s", tmp)==1) { // the one and only correct condition!!!
			if(strcmp(tmp, BORDER)) {
				continue;
			}
			else {
				fscanf(fp, "%d", (int*)p++ );
			} // end else
		} // end while fscanf

		fclose(fp);
		counter++;

	} // end while

	p = p0;


	qsort((void*)p, tot_bord, sizeof(Bord_Val), Cmp_Bord_Val);

	db=1;

	for(i=1; i<tot_bord; i++ ) {
		if( p[i] != p[i-1] ) db++;
	}


	counter=1;

	// cl->infile is a namespace here
   sprintf(file, "%s_%0*d.%s", cl->infile, SGM_SUFFIX_WIDTH, counter, EXT1);

   if( !Input_Segment(file, &s_tmp) )
      return EXIT_FAILURE;


   strcpy(s.source, s_tmp.source);
   strcpy(s.a.p, s_tmp.a.p);
   s.a.size   = s_tmp.a.size;
   s.function = s_tmp.function;
	s.num = db;

   Flush_Segment(&s_tmp);


   if(!Alloc_Borders(db, &(s.b)))
      return EXIT_FAILURE;

   if(!Read_Sequence(s.source, &seq))
      return EXIT_FAILURE;


	s.b[0].k=p[0];
	s.b[0].prob=0;

	j=1;

	for(i=1; i<tot_bord; i++) {
		if(p[i]!=p[i-1]) {
			s.b[j].k = p[i];
			s.b[j].prob = 0;
			Fill_Counts(seq, s.b[j-1].k, s.b[j-1].cnt, s.b[j].k, s.b[j].cnt);
			j++;
		} // end if
	} // end for

	// cl->infile is a namespace here
	strcpy(outfile, cl->infile);
	strcat(outfile, ".sgm");

	Print_Segment(outfile, &s);

	Flush_Segment(&s);
	return EXIT_SUCCESS;
} // end main()
//--------------------------------------------------------------------
static int Cmp_Bord_Val(void const *p1, void const *p2)
{
	const Bord_Val *b1 = p1, *b2 = p2;
   if(*b1<*b2)
      return -1;
   else if(*b1>*b2)
      return 1;

   return 0;

} // end Cmp_Bord_Val()
//--------------------------------------------------------------------
int  Process_Options(Comm_Line *c)
{

   if(strchr(c->options,'h') || strchr(c->options,'?'))
      //return NULL;
      return 0;

   return OK;

} // end Process_Options()
//--------------------------------------------------------------------
