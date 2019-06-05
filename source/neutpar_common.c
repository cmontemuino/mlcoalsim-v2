#include "neut_tests.h"

#include "zf_log.h"

int ispolnomhit(long int j,int init, int nsam, int totalsam, char **dna_matrix, long int *posit)
{
    int i,h,g;
    int a0,a1;
    
    if(j) if(posit[j-1] == posit[j]) {
		return -1;/*estem a la segona mutació o més*/
	}
    
	/*mhit, including outgroup*//**/
    h = g = a1 = 0;
	a0 = '0';
    for(i=0;i<totalsam;i++) {
        if((int)dna_matrix[i][j] != a0) {
            if(!a1) a1 = (int)dna_matrix[i][j];
            if((int)dna_matrix[i][j] != a1) {
				return -2;
			}
        }
    }
	/**/
    h = g = a1 = 0;
	a0 = '0';
    for(i=init;i<init+nsam;i++) {
        if((int)dna_matrix[i][j] == a0) h++;
        else {
            if(!a1) a1 = (int)dna_matrix[i][j];
            if((int)dna_matrix[i][j] == a1) g++;
            else {
				return(-2);/*en el cas de la primera mutacio hagin almenys tres variants*/
			}
        }
    }
    if(h==nsam || h == 0) {
		return 0; /*invariant intrapop*/
	}

    return(nsam-h); /*gives the frequency of the variant that is different of '0'*/
}
