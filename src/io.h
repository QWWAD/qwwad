/************************************/
/*  function reads character from   */
/*  standard input without echoing  */
/************************************/
char
getch()
{
char c;

system ("stty raw -echo");
c = getchar();
system ("stty cooked echo");
return (c);
}


char 
getche()
{
char c;

system("stty raw");
c=getchar();
system("stty cooked");

return c;
}

/*****************************************/
/*  function clears the terminal screen  */
/*****************************************/
void clrscr()
{
system("clear");
}


double   *read_file();

double *read_file(s,n,l)

/* This function reads a file into memory and returns the start
   address of this block of memory, a pointer to a double which can be 
   converted into a pointer to a structure by a cast, e.g.  
   start_address = (struct files *)read_v("v.r", &n, 3);
   where the arguments are s = filename, address of a string, e.g. "v.r",
                           n = number of lines, address of an integer
                           l = number of doubles per line, integer 
   The calloc function call is defined in malloc.h which can be included by
   #include <malloc.h>
   The memory area is initialized with zeros.
   The number of lines is changed and returned via its address.
   The allocated memory can be removed by
   free(start_address);                                            */

char	s[];		/* filename				   */
int	*n;		/* number of lines			   */
int	l;		/* number of elements per line		   */
{
int	i;		/* index				   */
int	j;		/* general integer			   */
FILE	*fp;		/* file pointer                            */
double	*start;		/* start address of allocated memory	   */
double	*p;		/* temporary pointer             	   */
 
  if((fp=fopen(s,"r"))==0)  {
    fprintf(stderr,"Error: Cannot open file %s!\n", s);
    exit(0);
  }

  *n=0;				/* count number of lines in file   */
  do  {
    for (i=0; i<l; i++)
      j = fscanf(fp,"%*le");
    (*n)++;
  }  while (j!=EOF);
  (*n)--;
  rewind(fp);

  start = (double *)calloc(*n,l*sizeof(double));		/* allocate memory   */
  if (start==0)  {
   fprintf(stderr,"Cannot allocate memory!\n");
   exit(0);
  }
  p = start;

  while(fscanf(fp,"%le", p)!=EOF)		/* write into memory   */
   p++;

  fclose(fp);
  return(start);
}
