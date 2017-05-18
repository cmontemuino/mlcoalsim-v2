#include <stdlib.h>
#include <math.h>

static long int idum = 0;

/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */

#define IA 16807
#define IM 2147483647
#define IQ 127773
#define IR 2836
#define NTAB 32
#define EPS (1.2E-07)
#define MAX(a,b) (a>b)?a:b
#define MIN(a,b) (a<b)?a:b

static long int iv[NTAB],iy=0;
static double NDIV = (double)1.0/((double)1.0+((double)IM-(double)1.0)/(double)NTAB);
static double RNMX = ((double)1.0-(double)EPS);
static double AM = ((double)1.0/(double)IM);

void init_seed1(long int seed)
{
    int j;
    
    idum = -seed;
    for(j=0;j<NTAB;j++) iv[j]=0;
    iy=0;
}

double ran1(void)
{
	long int j,k;
	double ran;

	if ((idum <= 0) || (iy == 0)) {
		idum = MAX(-idum,idum);
		for(j=NTAB+7;j>=0;j--) {
			k = (long int)((double)idum/(double)IQ);
			idum = IA*(idum-k*IQ)-IR*k;
			if(idum < 0) idum += IM;
			if(j < NTAB) iv[j] = idum;
		}
		iy = iv[0];
	}
	k = (long int)((double)idum/(double)IQ);
	idum = IA*(idum-k*IQ)-IR*k;
	if(idum<0) idum += IM;
	j = (long int)((double)iy*NDIV);
	iy = iv[j];
	iv[j] = idum;
	ran = MIN((double)AM*(double)iy,(double)RNMX);
	
	return (double)ran;
}

double randoms(void)
{
	int rand();
	return(rand()/(RAND_MAX+1.0));
}

